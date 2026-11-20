#include <spdlog/spdlog.h>
#include "representation/IfcGeometry.h"
#include "representation/geometry.h"
#include "nurbs.h"
#include "operations/geometryutils.h"
#include <tinynurbs/tinynurbs.h>
#include <CDT.h>
#include <numeric>
#include <string>

namespace webifc::geometry{

	// How far apart two neighbouring directions of the surface may be before another row of points goes in
	constexpr double TRIM_TURN_STEP { 0.35 };
	constexpr int TRIM_MAX_STEPS { 48 };
	constexpr int TRIM_SEED_GRID { 12 };

	void Nurbs::fill_geometry(){
		if (!_initialized) {
			return;
		}
		if (this->fill_geometry_trimmed()) return;
		this->fill_geometry_legacy();
	}

	Nurbs::uv_point_t Nurbs::wrapped(uv_point_t const& uv) const
	{
		const auto into = [](double x, double lo, double hi, bool wrap)
		{
			if (!wrap) return glm::clamp(x, lo, hi);
			if (x >= lo && x <= hi) return x;
			double r = std::fmod(x - lo, hi - lo);
			if (r < 0) r += hi - lo;
			return lo + r;
		};
		return { into(uv.x, range_knots_u.x, range_knots_u.y, wrapU), into(uv.y, range_knots_v.x, range_knots_v.y, wrapV) };
	}

	double Nurbs::project(glm::dvec3 const& pt, uv_point_t& uv) const
	{
		double best = std::numeric_limits<double>::max();
		uv_point_t at = uv, bestUV = uv;
		for (int it = 0; it < 24; it++)
		{
			const uv_point_t in = this->wrapped(at);
			const auto d = tinynurbs::surfaceDerivatives(*this->nurbs, 1, in.x, in.y);
			const glm::dvec3 r = d(0, 0) - pt, su = d(1, 0), sv = d(0, 1);
			const double dist = glm::length(r);
			if (dist < best) { best = dist; bestUV = at; }
			const double a = glm::dot(su, su), b = glm::dot(su, sv), c = glm::dot(sv, sv), e = -glm::dot(r, su), f = -glm::dot(r, sv);
			double du = 0, dv = 0;
			const double det = a * c - b * b;
			if (std::abs(det) > 1e-14 * std::max(1e-300, a * c)) { du = (e * c - b * f) / det; dv = (a * f - b * e) / det; }
			else if (a > c && a > 0) du = e / a; // one direction has collapsed (a pole): move along the other
			else if (c > 0) dv = f / c;
			const double periodU = range_knots_u.y - range_knots_u.x, periodV = range_knots_v.y - range_knots_v.x;
			const uv_point_t next {
				wrapU ? at.x + glm::clamp(du, -periodU / 4, periodU / 4) : glm::clamp(at.x + du, range_knots_u.x, range_knots_u.y),
				wrapV ? at.y + glm::clamp(dv, -periodV / 4, periodV / 4) : glm::clamp(at.y + dv, range_knots_v.x, range_knots_v.y) };
			const double moved = std::abs(next.x - at.x) / (range_knots_u.y - range_knots_u.x) + std::abs(next.y - at.y) / (range_knots_v.y - range_knots_v.x);
			at = next;
			if (moved < 1e-12) break;
		}
		const uv_point_t in = this->wrapped(at);
		const double last = glm::distance(tinynurbs::surfacePoint(*this->nurbs, in.x, in.y), pt);
		if (last < best) { best = last; bestUV = at; }
		uv = bestUV;
		return best;
	}

	bool Nurbs::fill_geometry_trimmed()
	{
		try
		{
			const double u0 = range_knots_u.x, u1 = range_knots_u.y, v0 = range_knots_v.x, v1 = range_knots_v.y;
			if (!(u1 > u0) || !(v1 > v0)) return false;
			const double close = 1e-4 / this->scaling; // a tenth of a millimetre
			const double lost = 5e-3 / this->scaling;

			// A coarse picture of the surface: where to start each search, and how much it turns
			const int G = TRIM_SEED_GRID;
			std::vector<glm::dvec3> grid((G + 1) * (G + 1));
			for (int i = 0; i <= G; i++)
				for (int j = 0; j <= G; j++)
					grid[i * (G + 1) + j] = tinynurbs::surfacePoint(*this->nurbs, u0 + (u1 - u0) * i / G, v0 + (v1 - v0) * j / G);
			const auto seed = [&](glm::dvec3 const& pt, double& interior)
			{
				double best = std::numeric_limits<double>::max(); int bi = 0, bj = 0;
				for (int i = 0; i <= G; i++)
					for (int j = 0; j <= G; j++)
					{
						const glm::dvec3 d = grid[i * (G + 1) + j] - pt; const double dd = glm::dot(d, d);
						if (dd < best) { best = dd; bi = i; bj = j; }
					}
				interior = std::min(std::min(bi, G - bi), std::min(bj, G - bj));
				return uv_point_t { u0 + (u1 - u0) * bi / G, v0 + (v1 - v0) * bj / G };
			};

			// Does the surface meet itself across a direction?
			{
				bool u = true, v = true; double extentU = 0, extentV = 0;
				for (int k = 0; k <= G; k++)
				{
					if (glm::distance(grid[k], grid[G * (G + 1) + k]) > close * 10) u = false;
					if (glm::distance(grid[k * (G + 1)], grid[k * (G + 1) + G]) > close * 10) v = false;
					extentU = std::max(extentU, glm::distance(grid[k], grid[(G / 2) * (G + 1) + k]));
					extentV = std::max(extentV, glm::distance(grid[k * (G + 1)], grid[k * (G + 1) + G / 2]));
				}
				this->wrapU = u && extentU > close * 100;
				this->wrapV = v && extentV > close * 100;
			}
			const double periodU = u1 - u0, periodV = v1 - v0;

			// Each bound followed round in the parameter plane, from where it is least in doubt
			struct Loop { std::vector<uv_point_t> uv; std::vector<glm::dvec3> at; int roundU = 0; int roundV = 0; };
			std::vector<Loop> loops;
			for (auto const& bound : this->bounds)
			{
				// a corner given by both edges that meet there is one point
				std::vector<glm::dvec3> pts;
				pts.reserve(bound.curve.points.size());
				for (auto const& pt : bound.curve.points) if (pts.empty() || glm::distance(pts.back(), pt) > close * 1e-2) pts.push_back(pt);
				size_t n = pts.size();
				while (n > 1 && glm::distance(pts[n - 1], pts[0]) <= close * 1e-2) n--;
				if (n < 3) { if (&bound == &this->bounds.front()) return false; continue; }
				size_t start = 0; double mostInterior = -1;
				std::vector<uv_point_t> seeds(n);
				for (size_t i = 0; i < n; i++)
				{
					double interior = 0; seeds[i] = seed(pts[i], interior);
					if (interior > mostInterior) { mostInterior = interior; start = i; }
				}
				const auto follow = [&](glm::dvec3 const& pt, uv_point_t const& previous, uv_point_t fresh, double& left)
				{
					uv_point_t uv = previous;
					left = this->project(pt, uv);
					if (left > close)
					{
						// a fresh start, on the turn of the surface nearest the last point
						if (this->wrapU) fresh.x += std::round((previous.x - fresh.x) / periodU) * periodU;
						if (this->wrapV) fresh.y += std::round((previous.y - fresh.y) / periodV) * periodV;
						const double freshLeft = this->project(pt, fresh);
						if (freshLeft < left) { left = freshLeft; uv = fresh; }
					}
					return uv;
				};
				Loop loop; loop.uv.resize(n); loop.at.resize(n);
				uv_point_t previous = seeds[start];
				for (size_t k = 0; k < n; k++)
				{
					const size_t i = (start + k) % n;
					double left = 0;
					const uv_point_t uv = follow(pts[i], previous, seeds[i], left);
					if (left > lost) return false;
					loop.uv[k] = uv; loop.at[k] = pts[i]; previous = uv;
				}
				// back at the start: on the same turn of the surface, or one further round?
				double left = 0;
				const uv_point_t back = follow(pts[start], previous, seeds[start], left);
				if (this->wrapU) loop.roundU = static_cast<int>(std::round((back.x - loop.uv[0].x) / periodU));
				if (this->wrapV) loop.roundV = static_cast<int>(std::round((back.y - loop.uv[0].y) / periodV));
				loops.push_back(std::move(loop));
			}
			if (loops.empty()) return false;

			std::vector<CDT::V2d<double>> vertices;
			std::vector<glm::dvec3> positions;
			std::vector<CDT::Edge> edges;
			const auto addVertex = [&](uv_point_t const& uv, glm::dvec3 const& at)
			{
				vertices.push_back(CDT::V2d<double>{ uv.x, uv.y });
				positions.push_back(at);
				return static_cast<CDT::VertInd>(vertices.size() - 1);
			};
			const auto centre = [](Loop const& loop)
			{
				uv_point_t c(0); for (auto const& uv : loop.uv) c += uv;
				return c / static_cast<double>(loop.uv.size());
			};
			std::vector<size_t> round, plain;
			for (size_t l = 0; l < loops.size(); l++) (loops[l].roundU != 0 || loops[l].roundV != 0 ? round : plain).push_back(l);
			uv_point_t home = centre(loops[plain.empty() ? round.front() : plain.front()]);
			if (!round.empty())
			{
				// A face that runs right round the surface: two bounds, each once round the same way. It is cut
				// open along a line from one bound to the other.
				if (round.size() != 2) return false;
				Loop& a = loops[round[0]]; Loop& b = loops[round[1]];
				const bool alongU = a.roundU != 0;
				if (alongU ? (std::abs(a.roundU) != 1 || std::abs(b.roundU) != 1 || a.roundV != 0 || b.roundV != 0)
				           : (std::abs(a.roundV) != 1 || std::abs(b.roundV) != 1 || a.roundU != 0 || b.roundU != 0)) return false;
				const double period = alongU ? periodU : periodV;
				const auto along = [&](uv_point_t const& uv) -> double { return alongU ? uv.x : uv.y; };
				const auto shift = [&](uv_point_t& uv, double by) { if (alongU) uv.x += by; else uv.y += by; };
				// both running the way the parameter grows
				const auto forwards = [&](Loop& loop)
				{
					if ((alongU ? loop.roundU : loop.roundV) > 0) return;
					std::reverse(loop.uv.begin(), loop.uv.end()); std::reverse(loop.at.begin(), loop.at.end());
				};
				forwards(a); forwards(b);
				// the second starts beside the start of the first
				{
					size_t nearest = 0; double least = std::numeric_limits<double>::max();
					for (size_t i = 0; i < b.uv.size(); i++)
					{
						double d = std::fmod(along(b.uv[i]) - along(a.uv[0]), period);
						if (d < 0) d += period;
						d = std::min(d, period - d);
						if (d < least) { least = d; nearest = i; }
					}
					for (size_t i = 0; i < nearest; i++) shift(b.uv[i], period);
					std::rotate(b.uv.begin(), b.uv.begin() + nearest, b.uv.end());
					std::rotate(b.at.begin(), b.at.begin() + nearest, b.at.end());
					const double by = std::round((along(a.uv[0]) - along(b.uv[0])) / period) * period;
					for (auto& uv : b.uv) shift(uv, by);
				}
				CDT::VertInd ends[2][2];
				int which = 0;
				for (Loop* loop : { &a, &b })
				{
					const auto first = addVertex(loop->uv[0], loop->at[0]);
					for (size_t i = 1; i < loop->uv.size(); i++) addVertex(loop->uv[i], loop->at[i]);
					uv_point_t again = loop->uv[0]; shift(again, period);
					const auto last = addVertex(again, loop->at[0]);
					for (CDT::VertInd v = first; v < last; v++) edges.emplace_back(v, v + 1);
					ends[which][0] = first; ends[which][1] = last; which++;
				}
				edges.emplace_back(ends[0][0], ends[1][0]);
				edges.emplace_back(ends[0][1], ends[1][1]);
				home = a.uv[0]; shift(home, period / 2);
			}
			for (size_t l : plain)
			{
				Loop& loop = loops[l];
				// every bound on the same turn of the surface
				const uv_point_t c = centre(loop);
				const double byU = this->wrapU ? std::round((home.x - c.x) / periodU) * periodU : 0, byV = this->wrapV ? std::round((home.y - c.y) / periodV) * periodV : 0;
				const auto first = static_cast<CDT::VertInd>(vertices.size());
				const auto n = static_cast<CDT::VertInd>(loop.uv.size());
				for (CDT::VertInd i = 0; i < n; i++)
				{
					addVertex({ loop.uv[i].x + byU, loop.uv[i].y + byV }, loop.at[i]);
					edges.emplace_back(first + i, first + (i + 1) % n);
				}
			}
			double minU = std::numeric_limits<double>::max(), maxU = std::numeric_limits<double>::lowest(), minV = minU, maxV = maxU;
			for (auto const& v : vertices) { minU = std::min(minU, v.x); maxU = std::max(maxU, v.x); minV = std::min(minV, v.y); maxV = std::max(maxV, v.y); }
			if (vertices.size() < 3 || !(maxU > minU) || !(maxV > minV)) return false;
			const size_t boundCount = vertices.size();

			// Rows of points inside the face: as many as the surface turns across it
			const auto turning = [&](bool alongU)
			{
				double most = 0;
				for (int line = 0; line <= G; line++)
				{
					double sum = 0; glm::dvec3 before(0); bool has = false;
					for (int k = 0; k < G; k++)
					{
						const double lo = alongU ? u0 + (u1 - u0) * k / G : v0 + (v1 - v0) * k / G, hi = alongU ? u0 + (u1 - u0) * (k + 1) / G : v0 + (v1 - v0) * (k + 1) / G;
						if (!(alongU ? this->wrapU : this->wrapV) && (hi < (alongU ? minU : minV) || lo > (alongU ? maxU : maxV))) continue;
						const glm::dvec3 step = alongU ? grid[(k + 1) * (G + 1) + line] - grid[k * (G + 1) + line] : grid[line * (G + 1) + k + 1] - grid[line * (G + 1) + k];
						const double len = glm::length(step); if (len < 1e-12) continue;
						const glm::dvec3 dir = step / len;
						if (has) sum += std::acos(glm::clamp(glm::dot(dir, before), -1.0, 1.0));
						before = dir; has = true;
					}
					most = std::max(most, sum);
				}
				// the coarse picture sees a turn only between its own steps: allow for one more
				const double span = std::min(1.0, alongU ? (maxU - minU) / (u1 - u0) : (maxV - minV) / (v1 - v0));
				return most * std::min(1.0, span * G / std::max(1.0, span * G - 1.0) );
			};
			const int stepsU = glm::clamp(static_cast<int>(std::ceil(turning(true) / TRIM_TURN_STEP)), 1, TRIM_MAX_STEPS);
			const int stepsV = glm::clamp(static_cast<int>(std::ceil(turning(false) / TRIM_TURN_STEP)), 1, TRIM_MAX_STEPS);
			const double stepU = (maxU - minU) / stepsU, stepV = (maxV - minV) / stepsV;
			for (int i = 1; i < stepsU; i++)
				for (int j = 1; j < stepsV; j++)
				{
					const double u = minU + stepU * i, v = minV + stepV * j;
					// not on top of the boundary: that only makes slivers
					bool clear = true;
					for (size_t e = 0; e < edges.size() && clear; e++)
					{
						const auto& a = vertices[edges[e].v1()]; const auto& b = vertices[edges[e].v2()];
						const double ax = (a.x - u) / stepU, ay = (a.y - v) / stepV, bx = (b.x - u) / stepU, by = (b.y - v) / stepV;
						const double dx = bx - ax, dy = by - ay, len2 = dx * dx + dy * dy;
						const double t = len2 > 0 ? glm::clamp(-(ax * dx + ay * dy) / len2, 0.0, 1.0) : 0.0;
						const double px = ax + dx * t, py = ay + dy * t;
						if (px * px + py * py < 0.09) clear = false;
					}
					if (!clear) continue;
					vertices.push_back(CDT::V2d<double>{ u, v });
					const uv_point_t in = this->wrapped({ u, v });
					positions.push_back(tinynurbs::surfacePoint(*this->nurbs, in.x, in.y));
				}

			// The parameter plane stretched to the surface's own distances: there the triangulation picks the
			// triangles that are well shaped on the surface
			{
				double alongU = 0, alongV = 0;
				for (int line = 0; line <= G; line++)
					for (int k = 0; k < G; k++)
					{
						alongU += glm::distance(grid[(k + 1) * (G + 1) + line], grid[k * (G + 1) + line]);
						alongV += glm::distance(grid[line * (G + 1) + k + 1], grid[line * (G + 1) + k]);
					}
				const double perU = alongU / (G + 1) / (u1 - u0), perV = alongV / (G + 1) / (v1 - v0);
				if (perU > 0 && perV > 0 && std::isfinite(perU) && std::isfinite(perV))
					for (auto& v : vertices) { v.x *= perU; v.y *= perV; }
			}

			// Bound points that land on one another in the parameter plane (a pole, a seam met twice) are one vertex
			const auto duplicates = CDT::RemoveDuplicatesAndRemapEdges(vertices, edges);
			if (!duplicates.duplicates.empty())
			{
				std::vector<glm::dvec3> kept(vertices.size());
				for (size_t i = 0; i < positions.size(); i++) kept[duplicates.mapping[i]] = positions[i];
				// the first of each group wins, as the remap keeps it
				for (size_t i = positions.size(); i-- > 0;) kept[duplicates.mapping[i]] = positions[i];
				positions = std::move(kept);
			}
			edges.erase(std::remove_if(edges.begin(), edges.end(), [](CDT::Edge const& e) { return e.v1() == e.v2(); }), edges.end());
			std::sort(edges.begin(), edges.end(), [](CDT::Edge const& a, CDT::Edge const& b) { return a.v1() != b.v1() ? a.v1() < b.v1() : a.v2() < b.v2(); });
			// an edge walked twice (a seam) bounds nothing
			std::vector<CDT::Edge> once;
			for (size_t i = 0; i < edges.size();)
			{
				size_t j = i; while (j < edges.size() && edges[j] == edges[i]) j++;
				if ((j - i) % 2 == 1) once.push_back(edges[i]);
				i = j;
			}
			if (once.size() < 3) return false;

			CDT::Triangulation<double> triangulator{ CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::NotAllowed, 0.0 };
			triangulator.insertVertices(vertices);
			triangulator.insertEdges(once);
			triangulator.eraseOuterTrianglesAndHoles();
			if (triangulator.triangles.empty()) return false;
			(void)boundCount;
			for (auto const& triangle : triangulator.triangles)
			{
				const auto i0 = triangle.vertices[0], i1 = triangle.vertices[1], i2 = triangle.vertices[2];
				const auto& a = vertices[i0]; const auto& b = vertices[i1]; const auto& c = vertices[i2];
				glm::dvec3 p0 = positions[i0], p1 = positions[i1], p2 = positions[i2];
				if (glm::length(glm::cross(p1 - p0, p2 - p0)) < 1e-18) continue;
				// Positive UV winding follows dS/du x dS/dv, regardless of bound order.
				if ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) < 0) std::swap(p1, p2);
				geometry.AddFace(p0, p1, p2);
			}
			return true;
		}
		catch (...)
		{
			return false;
		}
	}

	void Nurbs::fill_geometry_legacy(){
		auto uv_points {this->get_uv_points()};
		auto indices {get_triangulation_uv_points(uv_points)};

		// Subdivide resulting triangles to increase definition
		// r indicates the level of subdivision, currently 3 you can increase it to 5
		for (size_t r = 0; r < 3; r++)
		{
			auto num_indices{indices.size()};
			std::vector<uint32_t> newIndices;
			newIndices.reserve(num_indices / 3 * 12);
			Nurbs::uv_points_t newUVPoints;
			newUVPoints.reserve(num_indices / 3 * 6);

			for (size_t i = 0; i < num_indices; i += 3)
			{
				auto const& p0 = newUVPoints.emplace_back(std::move(uv_points[indices[i + 0]]));
				auto const& p1 = newUVPoints.emplace_back(std::move(uv_points[indices[i + 1]]));
				auto const& p2 = newUVPoints.emplace_back(std::move(uv_points[indices[i + 2]]));
				newUVPoints.emplace_back((p0.x + p1.x) / 2, (p0.y + p1.y) / 2);
				newUVPoints.emplace_back((p0.x + p2.x) / 2, (p0.y + p2.y) / 2);
				newUVPoints.emplace_back((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);

				int offset = newUVPoints.size() - 6;

				newIndices.push_back(offset + 0);
				newIndices.push_back(offset + 3);
				newIndices.push_back(offset + 4);

				newIndices.push_back(offset + 3);
				newIndices.push_back(offset + 5);
				newIndices.push_back(offset + 4);

				newIndices.push_back(offset + 3);
				newIndices.push_back(offset + 1);
				newIndices.push_back(offset + 5);

				newIndices.push_back(offset + 4);
				newIndices.push_back(offset + 5);
				newIndices.push_back(offset + 2);
			}

			uv_points = newUVPoints;
			indices = newIndices;
		}

		for (size_t i = 0; i < indices.size(); i += 3)
		{
			auto const& p0 {uv_points[indices[i + 0]]};
			auto const& p1 {uv_points[indices[i + 1]]};
			auto const& p2 {uv_points[indices[i + 2]]};
			auto pt00 {tinynurbs::surfacePoint(*this->nurbs, p0.x, p0.y)};
			auto pt01 {tinynurbs::surfacePoint(*this->nurbs, p1.x, p1.y)};
			auto pt10 {tinynurbs::surfacePoint(*this->nurbs, p2.x, p2.y)};
            // Positive UV winding follows dS/du x dS/dv, regardless of bound order.
            const auto u = p1 - p0, v = p2 - p0;
            if (u.x * v.y - u.y * v.x < 0) std::swap(pt01, pt10);
			geometry.AddFace(pt00, pt01, pt10);
		}
	}	

	Nurbs::Nurbs(IfcGeometry& geometry, std::vector<IfcBound3D>const & bounds, IfcSurface const& surface, double const scaling)
		: geometry{geometry}, 
			bounds{bounds}, 
			bspline_surface{surface.BSplineSurface},
			num_u{this->bspline_surface.ControlPoints.size()},
			num_v{this->bspline_surface.ControlPoints.front().size()},
			scaling{scaling} {
			this->init();
	}
	
	void Nurbs::init() {
		// Check that the control point grid has sufficient dimensions.
		// We need at least (degree + 1) control points in each direction.
		if (this->num_u < static_cast<size_t>(this->bspline_surface.UDegree) + 1) {
			spdlog::error("Insufficient control point rows: num_u = {} but UDegree = {} requires at least {} rows",
				this->num_u, this->bspline_surface.UDegree, this->bspline_surface.UDegree + 1);
			return; // Or throw an exception.
		}
		if (this->num_v < static_cast<size_t>(this->bspline_surface.VDegree) + 1) {
			spdlog::error("Insufficient control point columns: num_v = {} but VDegree = {} requires at least {} columns",
				this->num_v, this->bspline_surface.VDegree, this->bspline_surface.VDegree + 1);
			return; // Or throw an exception.
		}

		// Validate degrees.
		if (this->bspline_surface.UDegree < 0 || this->bspline_surface.VDegree < 0) {
			spdlog::error("Invalid degree values: UDegree={}, VDegree={}",
				this->bspline_surface.UDegree, this->bspline_surface.VDegree);
			return;
		}

		// Validate control points count.
		size_t expectedCount = this->num_u * this->num_v;
		auto controlPoints = this->get_control_points();
		if (controlPoints.size() != expectedCount) {
			spdlog::error("Control points count mismatch: expected {}, got {}", expectedCount, controlPoints.size());
			return;
		}
		auto weights = this->get_weights();
		if (weights.size() != expectedCount) {
			spdlog::error("Weights count mismatch: expected {}, got {}", expectedCount, weights.size());
			return;
		}

		// Create the NURBS surface.
		this->nurbs = std::make_shared<tinynurbs::RationalSurface3d>(
			static_cast<int>(this->bspline_surface.UDegree),
			static_cast<int>(this->bspline_surface.VDegree),
			this->get_knots(this->bspline_surface.UKnots, this->bspline_surface.UMultiplicity),
			this->get_knots(this->bspline_surface.VKnots, this->bspline_surface.VMultiplicity),
			tinynurbs::array2<glm::dvec3>{this->num_u, this->num_v, controlPoints},
			tinynurbs::array2<double>{this->num_u, this->num_v, weights}
		);

		// Check that the knot vectors are large enough.
		if (this->nurbs->knots_u.size() < static_cast<size_t>(this->nurbs->degree_u + 1)) {
			spdlog::error("Invalid knots_u: size {} is less than degree_u+1 ({})",
				this->nurbs->knots_u.size(), this->nurbs->degree_u + 1);
			return;
		}
		if (this->nurbs->knots_v.size() < static_cast<size_t>(this->nurbs->degree_v + 1)) {
			spdlog::error("Invalid knots_v: size {} is less than degree_v+1 ({})",
				this->nurbs->knots_v.size(), this->nurbs->degree_v + 1);
			return;
		}

		// Helper lambda to check if a knot vector is monotonic increasing.
		auto check_monotonic = [](const std::vector<double>& knots, const std::string& name) -> bool {
			for (size_t i = 1; i < knots.size(); i++) {
				if (knots[i] < knots[i - 1]) {
					spdlog::error("{} is not monotonic increasing at index {} ({} < {})", name, i, knots[i], knots[i - 1]);
					return false;
				}
			}
			return true;
			};

		if (!check_monotonic(this->nurbs->knots_u, "knots_u")) return;
		if (!check_monotonic(this->nurbs->knots_v, "knots_v")) return;

		// Ensure that we have enough knots to set the range.
		if (this->nurbs->knots_u.size() <= this->nurbs->degree_u ||
			this->nurbs->knots_u.size() <= this->nurbs->degree_u + 1) {
			spdlog::error("Not enough knots in knots_u to determine range, size={}, degree_u={}",
				this->nurbs->knots_u.size(), this->nurbs->degree_u);
			return;
		}
		if (this->nurbs->knots_v.size() <= this->nurbs->degree_v ||
			this->nurbs->knots_v.size() <= this->nurbs->degree_v + 1) {
			spdlog::error("Not enough knots in knots_v to determine range, size={}, degree_v={}",
				this->nurbs->knots_v.size(), this->nurbs->degree_v);
			return;
		}
		this->range_knots_u = {
			this->nurbs->knots_u[this->nurbs->degree_u],
			this->nurbs->knots_u[this->nurbs->knots_u.size() - this->nurbs->degree_u - 1]
		};
		this->range_knots_v = {
			this->nurbs->knots_v[this->nurbs->degree_v],
			this->nurbs->knots_v[this->nurbs->knots_v.size() - this->nurbs->degree_v - 1]
		};

		// Compute sample surface points.
		this->ptc = tinynurbs::surfacePoint(*this->nurbs, EPS_TINY, EPS_TINY);
		this->pth = tinynurbs::surfacePoint(*this->nurbs, 1.0, EPS_TINY);
		this->ptv = tinynurbs::surfacePoint(*this->nurbs, EPS_TINY, 1.0);

		// Compute distances for further use.
		this->dh = glm::distance(ptc, pth);
		this->dv = glm::distance(ptc, ptv);
		this->pr = (dh + 1) / (dv + 1);
		if (std::isnan(this->pr) || std::isinf(this->pr)) {
			this->pr = 1.0;
		}

		// Scale error tolerances.
		this->minError /= this->scaling;
		this->maxError /= this->scaling;
		_initialized = true;
	}
	std::vector<double> Nurbs::get_weights() const{
		std::vector<double> result(this->num_u * this->num_v);
		std::fill(result.begin(), result.end(), 1.0);
		// The loader fills WeightPoints, one row per row of control points
		auto const& rows = this->bspline_surface.WeightPoints;
		if (rows.size() != this->num_u) return result;
		for (size_t r{0}; r < rows.size(); ++r)
		{
			if (rows[r].size() != this->num_v) { std::fill(result.begin(), result.end(), 1.0); return result; }
			for (size_t i{0}; i < this->num_v; ++i) if (rows[r][i] > 0) result[r * this->num_v + i] = rows[r][i];
		}
		return result;
	}
	std::vector<double> Nurbs::get_knots(std::vector<double>const & bs_knots, std::vector<uint32_t> const & bs_mults) const{
		std::vector<double> result;
		auto knots_no_expanded {this->check_knots(bs_knots)}; 
		auto const num_srf_knots {std::accumulate(bs_mults.begin(), bs_mults.end(), 0.0)};
		result.reserve(num_srf_knots);
		for(size_t knot_i{0}; knot_i < bs_knots.size(); ++knot_i){
			auto const knot {knots_no_expanded[knot_i]};
			auto const knot_mult {bs_mults[knot_i]};
			for(size_t i{0}; i < knot_mult; ++i) result.push_back(knot);
		}
		return result;
	}
	std::vector<glm::dvec3> Nurbs::get_control_points() const{
		std::vector<glm::dvec3> result;
		size_t num_points{0};
		for(auto const& row : this->bspline_surface.ControlPoints) num_points += row.size();
		result.reserve(num_points);
		for(auto const& row : this->bspline_surface.ControlPoints)  std::copy(row.begin(), row.end(), std::back_inserter(result));
		return result;
	}
	Nurbs::uv_points_t Nurbs::get_uv_points() const{
		Nurbs::uv_points_t points;
		auto const& bound_points {this->bounds.front().curve.points};
		size_t num_points{bound_points.size()};
		points.resize(num_points);
		auto count{0};
		std::transform(bound_points.begin(), bound_points.end(), points.begin(), [&](auto const& point){
			++count;
			auto uv {this->inverse_evaluation(point)};
			return Nurbs::uv_point_t{uv.x, uv.y};				
		});
		std::sort(points.begin(), points.end(),[](auto const& left, auto const& right){
			  if (left[0] != right[0]) return left[0] < right[0];
        return left[1] < right[1];
		});
		auto last_it2 = std::unique(points.begin(), points.end(), [](auto const& a, auto const& b){
				double EPS{1E-5};
				return std::abs(a[0] - b[0]) < EPS && std::abs(a[1] - b[1]) < EPS;
		});
		points.erase(last_it2, points.end());
		std::sort(points.begin(), points.end(),[](auto const& left, auto const& right){
		  return left[1] < right[1];
		});
		return points;
	}
	auto Nurbs::get_approximation(glm::dvec3 const& pt, uv_point_t const& range_u, uv_point_t const& range_v) const{
		double fU{0.0};
		double fV{0.0};
    int const grid_size {10};
    auto min_distance = std::numeric_limits<double>::max();
		auto const portion_u {std::abs((range_u.y - range_u.x) / grid_size)};
		auto const portion_v {std::abs((range_v.y - range_v.x) / grid_size)};
		auto new_range_u {range_u};
		auto new_range_v {range_v};
		for (int i = 0; i < grid_size; ++i) {
				auto const step_u {portion_u * i};
				auto const u {range_u.x + step_u != 0.0 ? step_u : step_u + std::numeric_limits<double>::epsilon()};
        for (int j = 0; j < grid_size; ++j) {
					auto const step_v {portion_v * i};
					auto const v {range_v.x + step_v != 0.0 ? step_v : step_v + std::numeric_limits<double>::epsilon() };
					auto const pt_grid {tinynurbs::surfacePoint(*this->nurbs, u, v)};
					auto const dist {glm::distance(pt_grid, pt)};
					if (dist < min_distance) {
							min_distance = dist;
							fU = u;
							fV = v;
							new_range_u = {range_u.x + portion_u * i, range_u.x + portion_u * (i + 1) };
							new_range_v = {range_v.x + portion_v * i, range_v.x + portion_v * (i + 1) };
					}
        }
    }
		return std::tuple {min_distance, fU, fV, new_range_u, new_range_v};
	}	
	Nurbs::uv_point_t Nurbs::inverse_evaluation(glm::dvec3 const& pt) const
	{
		spdlog::debug("[BSplineInverseEvaluation({})]");
		return inverse_method(pt);
	}
	Nurbs::uv_point_t Nurbs::inverse_method(glm::dvec3 const& pt) const
	{
		spdlog::debug("[InverseMethod({})]");
		glm::highp_dvec3 pt00{};
		double fU {0.5};
		double fV {0.5};
		// auto [max_distance, fU, fV, new_range_u, new_range_v] {this->get_approximation(pt, this->range_knots_u, this->range_knots_v)};
		// if(max_distance <= maxError) return {fU, fV};
		// auto previous_distance{max_distance};
		// while(max_distance > maxError){
		// 	auto [next_max_distance, next_fU, next_fV, range_u, range_v] {this->get_approximation(pt, new_range_u, new_range_v)};
		// 	fU = next_fU;
		// 	fV = next_fV;
		// 	if(max_distance <= maxError) return {fU, fV};
		// 	if((previous_distance - next_max_distance) < 0.010) break;
		// 	previous_distance = next_max_distance;
		// 	new_range_u = range_u;
		// 	new_range_v = range_v;
		// 	max_distance = next_max_distance;
		// }
		
		size_t count{0};
		double divisor {100.0};
		auto max_distance = std::numeric_limits<double>::max();
		while (max_distance > maxError && divisor < 10000)
		{
			for (double r = 1; r < 5; r++)
			{
				int round = 0;
				auto mul_divisor {r * r * divisor};
				while (max_distance > minError && round < 3)
				{
					for (double i = 0; i < rotations; i++)
					{
						double rads = (i / rotations) * pi2;
						double incU = glm::sin(rads) / mul_divisor;
						double incV = glm::cos(rads) / mul_divisor;
						if (pr > 1) incV *= pr;
						else incU /= pr;
						while (true)
						{
							++count;
							// spdlog::debug("result: {} \t || {} \t || count: {} \t",fU, fV, count);
							double ffU = fU + incU;
							double ffV = fV + incV;
							if (ffU < range_knots_u.x)ffU = range_knots_u.y - (range_knots_u.x - ffU);
							else if (ffU > range_knots_u.y) ffU = range_knots_u.x + (ffU - range_knots_u.y);
							if (ffV < range_knots_v.x) ffV = range_knots_v.y - (range_knots_v.x - ffV);
							else if (ffV > range_knots_v.y) ffV = range_knots_v.x + (ffV - range_knots_v.y);
							pt00 = tinynurbs::surfacePoint(*this->nurbs, ffU, ffV);							
							auto const di {glm::distance(pt00, pt)};
							if (di < max_distance)
							{
								max_distance = di;
								fU = ffU;
								fV = ffV;
							}
							else
							{
								break;
							}
						}
					}
					round++;
				}
			}
			divisor *= 3;
		}
		return {fU, fV};
	}
	std::vector<uint32_t> Nurbs::get_triangulation_uv_points(Nurbs::uv_points_t const& uv_points) const{
		std::vector<uint32_t> result;
		if(uv_points.empty()) return result;
		auto const num_points {uv_points.size()};
		auto const num_edges {num_points};
		CDT::Triangulation<double> triangulator{};
		std::vector<CDT::V2d<double>> points;
		points.resize(num_points);
		std::transform(uv_points.begin(), uv_points.end(), points.begin(), [](auto const& uv_point){
			return CDT::V2d<double>{uv_point.x, uv_point.y};
		});

		try
		{
			triangulator.insertVertices(points);
			triangulator.eraseSuperTriangle();
			auto const num_indices{triangulator.triangles.size() * 3};
			result.reserve(triangulator.triangles.size() * 3);
			for(auto const& triangle : triangulator.triangles){
				auto const & vertice0_id{triangle.vertices[0]};
				auto const & vertice1_id{triangle.vertices[1]};
				auto const & vertice2_id{triangle.vertices[2]};
				auto const& vertice0 {glm::dvec2{uv_points[vertice0_id][0], uv_points[vertice0_id][1]}};
				auto const& vertice1 {glm::dvec2{uv_points[vertice1_id][0], uv_points[vertice1_id][1]}};
				auto const& vertice2 {glm::dvec2{uv_points[vertice2_id][0], uv_points[vertice2_id][1]}};
				auto const area{areaOfTriangle(vertice0, vertice1, vertice2)};
				constexpr double EPS {1E-2};
				constexpr double EPS2 {EPS*EPS};
				if(area < EPS2)
					continue;
				result.emplace_back(std::move(vertice0_id));
				result.emplace_back(std::move(vertice1_id));
				result.emplace_back(std::move(vertice2_id));
			}
		}
		catch(...){ return {};}
		return result;
	}
	std::vector<double> Nurbs::get_zscores(std::vector<double> const& knots) const{
		std::vector<double> result(knots.size());
		double mean = std::accumulate(knots.begin(), knots.end(), 0.0) / knots.size();
		double sq_sum = std::inner_product(knots.begin(), knots.end(), knots.begin(), 0.0);
		double stdev = std::sqrt(sq_sum / knots.size() - mean * mean);
		for (size_t i = 0; i < knots.size(); ++i) {
				result[i] = (knots[i] - mean) / stdev;
		}
		return result;
	}
	std::vector<double> Nurbs::check_knots(std::vector<double> const& knots) const{
		std::vector<double> result(knots.size());
		auto const num_knots {knots.size()};
		if(num_knots == 2){
			result[0] = 0.0;
			result[1] = 1.0;
			return result;
		}
		auto threshold {3.0};
		auto zscores {get_zscores(knots)};
		for(size_t i{0}; i < num_knots; ++i){
			if(std::abs(zscores[i]) > threshold){
				if(i == 0)									result[i] = knots[i+1];
				else if (i == num_knots -1)	result[i] = knots[i-1];
				else 												result[i] = (knots[i-1]+knots[i+1]) / 2.0;
			}
			else result[i] = knots[i]; 
		}
		return result;
	}
}