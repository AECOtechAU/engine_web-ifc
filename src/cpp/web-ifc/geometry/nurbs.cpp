#include <spdlog/spdlog.h>
#include "representation/IfcGeometry.h"
#include "representation/geometry.h"
#include "nurbs.h"
#include "operations/geometryutils.h"
#include <tinynurbs/tinynurbs.h>
#include <CDT.h>
#include <numeric>
#include <unordered_set>
#include <map>
#include <array>
#include <string>

namespace webifc::geometry{

	// How far apart two neighbouring directions of the surface may be before another row of points goes in
	// Free-form faces: rows of points every 12 degrees of turn, edges kept within a chord across 14
	constexpr double TRIM_FREEFORM_TURN_STEP { 0.21 };
	constexpr double TRIM_FREEFORM_STRAY { 0.03 };
	constexpr int TRIM_MAX_STEPS { 48 };
	constexpr int TRIM_SEED_GRID { 12 };
	constexpr int TRIM_REFINE_ROUNDS { 6 };
	// Where a surface turns evenly, an edge shorter than this many allowed turns (squared) is not measured:
	// its outline is sampled at the loader's circle steps, and edges across one such step are as good as the outline
	constexpr double TRIM_EVEN_SKIP { 4.0 };

	void Nurbs::fill_geometry(){
		if (!_initialized) {
			return;
		}
		if (this->fill_geometry_trimmed()) return;
		this->fill_geometry_legacy();
	}

	Nurbs::uv_point_t Nurbs::wrapped(uv_point_t const& uv) const
	{
		TrimSurface range;
		range.u0 = range_knots_u.x; range.u1 = range_knots_u.y; range.v0 = range_knots_v.x; range.v1 = range_knots_v.y;
		range.wrapU = wrapU; range.wrapV = wrapV;
		return range.wrapped(uv);
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

	namespace
	{
		struct Ball { std::vector<glm::dvec3> corners; std::vector<std::array<uint32_t, 3>> faces; };

		// The icosahedron with every edge cut into so many pieces and the new corners pushed out to the sphere
		Ball const& GeodesicBall(int cuts)
		{
			static std::map<int, Ball> made;
			auto found = made.find(cuts);
			if (found != made.end()) return found->second;
			const double t = (1.0 + std::sqrt(5.0)) / 2.0;
			const glm::dvec3 base[12] = { {-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1} };
			const int sides[20][3] = { {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
			                           {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1} };
			Ball ball;
			std::map<std::array<int64_t, 3>, uint32_t> known;
			const auto corner = [&](glm::dvec3 const& at)
			{
				const glm::dvec3 on = glm::normalize(at);
				const std::array<int64_t, 3> key { std::llround(on.x * 1e9), std::llround(on.y * 1e9), std::llround(on.z * 1e9) };
				auto it = known.find(key);
				if (it != known.end()) return it->second;
				const uint32_t index = static_cast<uint32_t>(ball.corners.size());
				ball.corners.push_back(on); known.emplace(key, index);
				return index;
			};
			for (auto const& side : sides)
			{
				const glm::dvec3 a = glm::normalize(base[side[0]]), b = glm::normalize(base[side[1]]), c = glm::normalize(base[side[2]]);
				const auto at = [&](int i, int j) { return corner(a * static_cast<double>(cuts - i - j) + b * static_cast<double>(i) + c * static_cast<double>(j)); };
				for (int i = 0; i < cuts; i++)
					for (int j = 0; j < cuts - i; j++)
					{
						ball.faces.push_back({ at(i, j), at(i + 1, j), at(i, j + 1) });
						if (i + j + 1 < cuts) ball.faces.push_back({ at(i + 1, j), at(i + 1, j + 1), at(i, j + 1) });
					}
			}
			return made.emplace(cuts, std::move(ball)).first->second;
		}
	}

	int BallCuts(double acrossMetres)
	{
		if (acrossMetres < 0.04) return 1;
		if (acrossMetres < 0.1) return 2;
		if (acrossMetres < 0.3) return 3;
		int cuts = 4;
		// a face's middle sits about 0.29 / cuts^2 of the radius inside the sphere
		while (cuts < 16 && 0.29 * (acrossMetres / 2) / (cuts * cuts) > 0.010) cuts++;
		return cuts;
	}

	void AddWholeBall(IfcGeometry& geometry, glm::dvec3 const& centre, double radius, double scaling)
	{
		Ball const& ball = GeodesicBall(BallCuts(2 * radius * scaling));
		for (auto const& face : ball.faces)
			geometry.AddFace(centre + ball.corners[face[0]] * radius, centre + ball.corners[face[1]] * radius, centre + ball.corners[face[2]] * radius);
	}

	bool TriangulateSphereFace(IfcGeometry& geometry, std::vector<IfcBound3D> const& bounds, glm::dvec3 const& centre, double radius, double scaling, int outerTurn, int looks)
	{
		try
		{
			if (!(radius > 0)) return false;
			const double close = 1e-4 / scaling;
			std::vector<std::vector<glm::dvec3>> loops;
			for (auto const& bound : bounds)
			{
				std::vector<glm::dvec3> pts;
				for (auto const& pt : bound.curve.points) if (pts.empty() || glm::distance(pts.back(), pt) > close * 1e-2) pts.push_back(pt);
				while (pts.size() > 1 && glm::distance(pts.back(), pts.front()) <= close * 1e-2) pts.pop_back();
				if (pts.size() < 3) { if (&bound == &bounds.front()) return false; continue; }
				for (auto const& pt : pts) if (std::abs(glm::distance(pt, centre) - radius) > std::max(5e-3 / scaling, radius * 0.02)) return false;
				loops.push_back(std::move(pts));
			}
			if (loops.empty()) return false;

			// Which side of its outline the face is on: the outline runs anticlockwise about the face's normal,
			// which looks outwards when the face looks the way the sphere does
			glm::dvec3 turn(0), mean(0);
			auto const& outer = loops.front();
			for (size_t i = 0; i < outer.size(); i++)
			{
				turn += glm::cross(outer[i] - centre, outer[(i + 1) % outer.size()] - centre);
				mean += outer[i] - centre;
			}
			mean /= static_cast<double>(outer.size());
			glm::dvec3 towards;
			if (outerTurn != 0 && glm::length(turn) > radius * radius * 1e-6) towards = glm::normalize(turn) * static_cast<double>(outerTurn * looks);
			else if (glm::length(mean) > radius * 0.05) towards = glm::normalize(mean);
			else return false;

			// The sphere laid flat from the point opposite the face: circles stay circles, so the triangles the
			// plane picks are the ones the sphere would
			glm::dvec3 e1 = glm::cross(towards, std::abs(towards.x) < 0.9 ? glm::dvec3(1, 0, 0) : glm::dvec3(0, 1, 0));
			e1 = glm::normalize(e1);
			const glm::dvec3 e2 = glm::cross(towards, e1);
			bool behind = false;
			const auto flat = [&](glm::dvec3 const& at)
			{
				const glm::dvec3 x = (at - centre) / radius;
				const double d = 1 + glm::dot(x, towards);
				if (d < 1e-3) { behind = true; return CDT::V2d<double>{ 0, 0 }; }
				return CDT::V2d<double>{ glm::dot(x, e1) / d, glm::dot(x, e2) / d };
			};
			std::vector<CDT::V2d<double>> vertices;
			std::vector<glm::dvec3> positions;
			std::vector<CDT::Edge> edges;
			for (auto const& loop : loops)
			{
				const auto first = static_cast<CDT::VertInd>(vertices.size());
				const auto n = static_cast<CDT::VertInd>(loop.size());
				for (CDT::VertInd i = 0; i < n; i++)
				{
					vertices.push_back(flat(loop[i])); positions.push_back(loop[i]);
					edges.emplace_back(first + i, first + (i + 1) % n);
				}
			}
			if (behind) return false;
			const size_t outline = vertices.size();

			// The ball's own corners that fall inside the face, clear of its outline
			const int cuts = BallCuts(2 * radius * scaling);
			Ball const& ball = GeodesicBall(cuts);
			const double clear = radius * 1.107 / cuts * 0.4;
			for (auto const& unit : ball.corners)
			{
				if (glm::dot(unit, towards) < -0.9) continue;
				const glm::dvec3 at = centre + unit * radius;
				const CDT::V2d<double> q = flat(at);
				bool inside = false, tooNear = false;
				for (size_t e = 0; e < outline && !tooNear; e++)
				{
					const auto& a = vertices[edges[e].v1()]; const auto& b = vertices[edges[e].v2()];
					if ((a.y > q.y) != (b.y > q.y) && q.x < (b.x - a.x) * (q.y - a.y) / (b.y - a.y) + a.x) inside = !inside;
					const glm::dvec3 pa = positions[edges[e].v1()], ab = positions[edges[e].v2()] - pa;
					const double len2 = glm::dot(ab, ab);
					const double t = len2 > 0 ? glm::clamp(glm::dot(at - pa, ab) / len2, 0.0, 1.0) : 0.0;
					if (glm::distance(at, pa + ab * t) < clear) tooNear = true;
				}
				if (!inside || tooNear) continue;
				vertices.push_back(q); positions.push_back(at);
			}
			behind = false;

			const auto duplicates = CDT::RemoveDuplicatesAndRemapEdges(vertices, edges);
			if (!duplicates.duplicates.empty())
			{
				std::vector<glm::dvec3> kept(vertices.size());
				for (size_t i = positions.size(); i-- > 0;) kept[duplicates.mapping[i]] = positions[i];
				positions = std::move(kept);
			}
			edges.erase(std::remove_if(edges.begin(), edges.end(), [](CDT::Edge const& e) { return e.v1() == e.v2(); }), edges.end());
			std::sort(edges.begin(), edges.end(), [](CDT::Edge const& a, CDT::Edge const& b) { return a.v1() != b.v1() ? a.v1() < b.v1() : a.v2() < b.v2(); });
			std::vector<CDT::Edge> once;
			for (size_t i = 0; i < edges.size();)
			{
				size_t j = i; while (j < edges.size() && edges[j] == edges[i]) j++;
				if ((j - i) % 2 == 1) once.push_back(edges[i]);
				i = j;
			}
			if (once.size() < 3) return false;
			CDT::Triangulation<double> whole{ CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::NotAllowed, 0.0 };
			whole.insertVertices(vertices);
			whole.insertEdges(once);
			whole.eraseOuterTrianglesAndHoles();
			if (whole.triangles.empty()) return false;
			for (auto const& triangle : whole.triangles)
			{
				glm::dvec3 p0 = positions[triangle.vertices[0]], p1 = positions[triangle.vertices[1]], p2 = positions[triangle.vertices[2]];
				const glm::dvec3 normal = glm::cross(p1 - p0, p2 - p0);
				if (glm::length(normal) < 1e-18) continue;
				// looking the way the surface does: out of the ball, or into it
				if ((glm::dot(normal, (p0 + p1 + p2) / 3.0 - centre) < 0) != (looks < 0)) std::swap(p1, p2);
				geometry.AddFace(p0, p1, p2);
			}
			return true;
		}
		catch (...)
		{
			return false;
		}
	}

	glm::dvec2 TrimSurface::wrapped(glm::dvec2 const& uv) const
	{
		const auto into = [](double x, double lo, double hi, bool wrap)
		{
			if (!wrap) return glm::clamp(x, lo, hi);
			if (x >= lo && x <= hi) return x;
			double r = std::fmod(x - lo, hi - lo);
			if (r < 0) r += hi - lo;
			return lo + r;
		};
		return { into(uv.x, u0, u1, wrapU), into(uv.y, v0, v1, wrapV) };
	}

	bool TriangulateTrimmedSurface(IfcGeometry& geometry, std::vector<IfcBound3D> const& bounds, TrimSurface& s, double scaling)
	{
		using uv_point_t = glm::dvec2;
		try
		{
			const double u0 = s.u0, u1 = s.u1, v0 = s.v0, v1 = s.v1;
			if (!(u1 > u0) || !(v1 > v0)) return false;
			const double close = 1e-4 / scaling; // a tenth of a millimetre
			const double lost = 5e-3 / scaling;

			// A coarse picture of the surface: where to start each search, and how much it turns
			const int G = s.direct ? TRIM_SEED_GRID / 2 : TRIM_SEED_GRID; // a surface worked out directly needs only its turning and lengths from this
			std::vector<glm::dvec3> grid((G + 1) * (G + 1));
			for (int i = 0; i <= G && !s.plain; i++)
				for (int j = 0; j <= G; j++)
					grid[i * (G + 1) + j] = s.point(u0 + (u1 - u0) * i / G, v0 + (v1 - v0) * j / G);
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
			if (s.findWraps)
			{
				bool u = true, v = true; double extentU = 0, extentV = 0;
				for (int k = 0; k <= G; k++)
				{
					if (glm::distance(grid[k], grid[G * (G + 1) + k]) > close * 10) u = false;
					if (glm::distance(grid[k * (G + 1)], grid[k * (G + 1) + G]) > close * 10) v = false;
					extentU = std::max(extentU, glm::distance(grid[k], grid[(G / 2) * (G + 1) + k]));
					extentV = std::max(extentV, glm::distance(grid[k * (G + 1)], grid[k * (G + 1) + G / 2]));
				}
				s.wrapU = u && extentU > close * 100;
				s.wrapV = v && extentV > close * 100;
			}
			const double periodU = u1 - u0, periodV = v1 - v0;

			// Where a surface closed round its first parameter comes to a point at an end of the second (a pole)
			bool poleLow = false, poleHigh = false;
			const glm::dvec3 tipLow = grid[0], tipHigh = grid[G];
			if (s.wrapU && !s.plain)
			{
				poleLow = poleHigh = true;
				for (int k = 1; k <= G; k++)
				{
					if (glm::distance(grid[k * (G + 1)], tipLow) > close * 10) poleLow = false;
					if (glm::distance(grid[k * (G + 1) + G], tipHigh) > close * 10) poleHigh = false;
				}
			}

			bool reachesLow = false, reachesHigh = false;
			// Each bound followed round in the parameter plane, from where it is least in doubt
			struct Loop { std::vector<uv_point_t> uv; std::vector<glm::dvec3> at; int roundU = 0; int roundV = 0; };
			std::vector<Loop> loops;
			for (auto const& bound : bounds)
			{
				// a corner given by both edges that meet there is one point
				std::vector<glm::dvec3> pts;
				pts.reserve(bound.curve.points.size());
				for (auto const& pt : bound.curve.points) if (pts.empty() || glm::distance(pts.back(), pt) > close * 1e-2) pts.push_back(pt);
				size_t n = pts.size();
				while (n > 1 && glm::distance(pts[n - 1], pts[0]) <= close * 1e-2) n--;
				if (n < 3) { if (&bound == &bounds.front()) return false; continue; }
				size_t start = 0; double mostInterior = -1;
				std::vector<uv_point_t> seeds(n);
				for (size_t i = 0; i < n && !s.direct; i++)
				{
					double interior = 0; seeds[i] = seed(pts[i], interior);
					if (interior > mostInterior) { mostInterior = interior; start = i; }
				}
				if (s.direct) seeds[0] = { (u0 + u1) / 2, (v0 + v1) / 2 };
				const auto follow = [&](glm::dvec3 const& pt, uv_point_t const& previous, uv_point_t fresh, double& left)
				{
					uv_point_t uv = previous;
					left = s.project(pt, uv);
					if (left > close && !s.direct)
					{
						// a fresh start, on the turn of the surface nearest the last point
						if (s.wrapU) fresh.x += std::round((previous.x - fresh.x) / periodU) * periodU;
						if (s.wrapV) fresh.y += std::round((previous.y - fresh.y) / periodV) * periodV;
						const double freshLeft = s.project(pt, fresh);
						if (freshLeft < left) { left = freshLeft; uv = fresh; }
					}
					return uv;
				};
				// A bound through a pole of the surface, or two (the outline of a slice of a ball): at a pole the turn
				// about the axis is free, so each run between poles is followed on its own, and the runs are set
				// side by side the way the bound's direction asks
				if (s.wrapU && (poleLow || poleHigh))
				{
					std::vector<int> atPole(n, 0);
					for (size_t i = 0; i < n; i++)
					{
						if (poleLow && glm::distance(pts[i], tipLow) < close * 10) atPole[i] = 1;
						else if (poleHigh && glm::distance(pts[i], tipHigh) < close * 10) atPole[i] = 2;
					}
					size_t poleRuns = 0, first = n;
					for (size_t i = 0; i < n; i++)
						if (atPole[i] && !atPole[(i + n - 1) % n]) { poleRuns++; if (first == n) first = i; }
					if (poleRuns >= 1)
					{
						struct Run { std::vector<uv_point_t> uv; std::vector<glm::dvec3> at; int before = 0, after = 0; };
						std::vector<Run> runs;
						int lastPole = 0;
						for (size_t k = 0; k < n;)
						{
							const size_t i = (first + k) % n;
							if (atPole[i]) { lastPole = atPole[i]; k++; continue; }
							Run run; run.before = lastPole;
							uv_point_t previous = s.direct ? (runs.empty() ? seeds[0] : runs.back().uv.back()) : seeds[i];
							for (; k < n && !atPole[(first + k) % n]; k++)
							{
								const size_t j = (first + k) % n;
								double left = 0;
								const uv_point_t uv = follow(pts[j], previous, s.direct ? previous : seeds[j], left);
								if (left > lost) return false;
								run.uv.push_back(uv); run.at.push_back(pts[j]); previous = uv;
							}
							run.after = atPole[(first + k) % n];
							runs.push_back(std::move(run));
						}
						if (runs.empty() || runs.size() > 4) return false;
						for (auto const& run : runs) { if (run.before == 1 || run.after == 1) reachesLow = true; if (run.before == 2 || run.after == 2) reachesHigh = true; }
						const int want = s.outerTurn * (&bound == &bounds.front() ? 1 : -1);
						const auto poleAt = [&](int which) { return which == 1 ? v0 : v1; };
						const auto build = [&](std::vector<int> const& turns, Loop& made)
						{
							made.uv.clear(); made.at.clear();
							for (size_t r = 0; r < runs.size(); r++)
							{
								const double by = turns[r] * periodU;
								// between a pole and a run that does not start beside it (the file left that edge out)
								// the outline follows the surface straight to the pole
								const auto leg = [&](uv_point_t const& from, uv_point_t const& to, glm::dvec3 const& atFrom, glm::dvec3 const& atTo)
								{
									if (glm::distance(atFrom, atTo) < close * 20) return;
									// as many pieces as the edge next door gets: an arc of this size and sweep
									const uv_point_t half = s.wrapped({ from.x, (from.y + to.y) / 2 });
									const glm::dvec3 middle = s.point(half.x, half.y);
									const double chord = glm::distance(atFrom, atTo);
									const double rise = glm::distance(middle, (atFrom + atTo) * 0.5);
									int pieces = 1;
									if (rise > chord * 1e-4)
									{
										const double radius = (chord * chord / 4 + rise * rise) / (2 * rise);
										pieces = glm::clamp(RoundSegments(radius * scaling, 2 * std::asin(glm::clamp(chord / (2 * radius), 0.0, 1.0))), 1, 64);
									}
									for (int piece = 1; piece < pieces; piece++)
									{
										const uv_point_t q { from.x, from.y + (to.y - from.y) * piece / pieces };
										const uv_point_t inside = s.wrapped(q);
										made.uv.push_back(q); made.at.push_back(s.point(inside.x, inside.y));
									}
								};
								const glm::dvec3 tipBefore = runs[r].before == 1 ? tipLow : tipHigh, tipAfter = runs[r].after == 1 ? tipLow : tipHigh;
								const uv_point_t poleBefore { runs[r].uv.front().x + by, poleAt(runs[r].before) }, poleAfter { runs[r].uv.back().x + by, poleAt(runs[r].after) };
								const uv_point_t first { runs[r].uv.front().x + by, runs[r].uv.front().y }, lastOf { runs[r].uv.back().x + by, runs[r].uv.back().y };
								made.uv.push_back(poleBefore); made.at.push_back(tipBefore);
								leg(poleBefore, first, tipBefore, runs[r].at.front());
								for (size_t q = 0; q < runs[r].uv.size(); q++) { made.uv.push_back({ runs[r].uv[q].x + by, runs[r].uv[q].y }); made.at.push_back(runs[r].at[q]); }
								leg(lastOf, poleAfter, runs[r].at.back(), tipAfter);
								made.uv.push_back(poleAfter); made.at.push_back(tipAfter);
							}
							double area = 0, along = 0;
							for (size_t q = 0; q < made.uv.size(); q++)
							{
								const uv_point_t& a = made.uv[q]; const uv_point_t& b = made.uv[(q + 1) % made.uv.size()];
								area += a.x * b.y - b.x * a.y;
								if (glm::distance(made.at[q], made.at[(q + 1) % made.uv.size()]) < close * 20) along += std::abs(a.x - b.x);
							}
							return std::pair<double, double>(area / 2, along);
						};
						// every way of setting the runs within a turn of one another
						std::vector<int> turns(runs.size(), 0), bestTurns;
						double bestAlong = std::numeric_limits<double>::max(), otherAlong = bestAlong; int bestSign = 0, otherSign = 0;
						Loop made;
						size_t combos = 1; for (size_t r = 1; r < runs.size(); r++) combos *= 3;
						for (size_t c = 0; c < combos; c++)
						{
							size_t rest = c;
							for (size_t r = 1; r < runs.size(); r++) { turns[r] = static_cast<int>(rest % 3) - 1; rest /= 3; }
							const auto [area, along] = build(turns, made);
							const int sign = area > 0 ? 1 : (area < 0 ? -1 : 0);
							if (sign == 0 || (want != 0 && runs.size() > 1 && sign != want)) continue;
							if (along < bestAlong - periodU * 1e-6) { otherAlong = bestAlong; otherSign = bestSign; bestAlong = along; bestSign = sign; bestTurns = turns; }
							else if (along < otherAlong) { otherAlong = along; otherSign = sign; }
						}
						if (bestTurns.empty()) return false;
						// without a direction, two equally short layouts on opposite sides cannot be told apart
						if (want == 0 && otherSign != 0 && otherSign != bestSign && otherAlong - bestAlong < periodU * 1e-3) return false;
						Loop loop;
						build(bestTurns, loop);
						loops.push_back(std::move(loop));
						continue;
					}
				}
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
				if (s.wrapU) loop.roundU = static_cast<int>(std::round((back.x - loop.uv[0].x) / periodU));
				if (s.wrapV) loop.roundV = static_cast<int>(std::round((back.y - loop.uv[0].y) / periodV));
				// a jump across both parameters whose corner is a pole of the surface goes through that corner
				for (size_t k = 0; k < loop.uv.size(); k++)
				{
					const size_t next = (k + 1) % loop.uv.size();
					uv_point_t a = loop.uv[k], b = loop.uv[next];
					if (next == 0) { b.x += loop.roundU * periodU; b.y += loop.roundV * periodV; }
					if (std::abs(a.x - b.x) < periodU * 0.02 || std::abs(a.y - b.y) < periodV * 0.02) continue;
					if (glm::distance(loop.at[k], loop.at[next]) < close * 100) continue;
					for (int which = 0; which < 2; which++)
					{
						const uv_point_t corner = which == 0 ? uv_point_t{ a.x, b.y } : uv_point_t{ b.x, a.y };
						const uv_point_t in = s.wrapped(corner);
						const glm::dvec3 at = s.point(in.x, in.y);
						const bool isA = glm::distance(at, loop.at[k]) < close * 10, isB = glm::distance(at, loop.at[next]) < close * 10;
						if (!isA && !isB) continue;
						// the way from a to the corner and on to b, with points where the surface turns
						std::vector<uv_point_t> pathUV; std::vector<glm::dvec3> pathAt;
						const glm::dvec3 atCorner = isA ? loop.at[k] : loop.at[next];
						const auto leg = [&](uv_point_t const& from, uv_point_t const& to, glm::dvec3 const& atFrom, glm::dvec3 const& atTo)
						{
							if (glm::distance(atFrom, atTo) < close * 10) return;
							const bool alongU = std::abs(from.x - to.x) > std::abs(from.y - to.y);
							const double start = alongU ? from.x : from.y, end = alongU ? to.x : to.y;
							std::vector<double> stops;
							if (s.creased == (alongU ? 0 : 1))
							{
								if (end > start) for (double c = std::floor(start) + 1; c < end - 1e-9; c += 1) stops.push_back(c);
								else for (double c = std::ceil(start) - 1; c > end + 1e-9; c -= 1) stops.push_back(c);
							}
							else
							{
								const int pieces = std::max(1, static_cast<int>(std::ceil(std::abs(end - start) / ((alongU ? periodU : periodV) / 16))));
								for (int piece = 1; piece < pieces; piece++) stops.push_back(start + (end - start) * piece / pieces);
							}
							for (double stop : stops)
							{
								const uv_point_t q = alongU ? uv_point_t{ stop, from.y } : uv_point_t{ from.x, stop };
								const uv_point_t inside = s.wrapped(q);
								pathUV.push_back(q); pathAt.push_back(s.point(inside.x, inside.y));
							}
						};
						const glm::dvec3 atA = loop.at[k], atB = loop.at[next];
						leg(a, corner, atA, atCorner);
						pathUV.push_back(corner); pathAt.push_back(atCorner);
						leg(corner, b, atCorner, atB);
						loop.uv.insert(loop.uv.begin() + k + 1, pathUV.begin(), pathUV.end());
						loop.at.insert(loop.at.begin() + k + 1, pathAt.begin(), pathAt.end());
						k += pathUV.size() - 1;
						k++;
						break;
					}
				}
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
				if (round.size() > 2) return false;
				Loop pole;
				if (round.size() == 1)
				{
					// One bound only: the face runs up to a point where the surface closes (a pole), which
					// stands in for the second bound
					Loop const& only = loops[round[0]];
					const bool u = only.roundU != 0;
					const auto closes = [&](bool low)
					{
						const int at = low ? 0 : G;
						for (int k = 1; k <= G; k++)
							if (glm::distance(u ? grid[k * (G + 1) + at] : grid[at * (G + 1) + k], u ? grid[at] : grid[at * (G + 1)]) > close * 10) return false;
						return true;
					};
					const bool low = closes(true), high = closes(false);
					if (low == high) return false;
					if (u) { if (low) reachesLow = true; else reachesHigh = true; }
					const double side = u ? (low ? v0 : v1) : (low ? u0 : u1);
					const glm::dvec3 tip = u ? grid[low ? 0 : G] : grid[(low ? 0 : G) * (G + 1)];
					pole.roundU = only.roundU; pole.roundV = only.roundV;
					for (int k = 0; k < 2; k++)
					{
						const double along = (u ? only.uv[0].x : only.uv[0].y) + (u ? periodU : periodV) * 0.5 * k * (((u ? only.roundU : only.roundV) > 0) ? 1 : -1);
						pole.uv.push_back(u ? uv_point_t{ along, side } : uv_point_t{ side, along });
						pole.at.push_back(tip);
					}
				}
				Loop& a = loops[round[0]]; Loop& b = round.size() == 2 ? loops[round[1]] : pole;
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
				const double byU = s.wrapU ? std::round((home.x - c.x) / periodU) * periodU : 0, byV = s.wrapV ? std::round((home.y - c.y) / periodV) * periodV : 0;
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
			const auto turning = [&](bool alongU, bool whole = false)
			{
				if (s.plain)
				{
					const double all = alongU ? s.turnU : s.turnV;
					return whole ? all : all * std::min(1.0, alongU ? (maxU - minU) / (u1 - u0) : (maxV - minV) / (v1 - v0));
				}
				double most = 0;
				for (int line = 0; line <= G; line++)
				{
					double sum = 0; glm::dvec3 before(0); bool has = false;
					for (int k = 0; k < G; k++)
					{
						const double lo = alongU ? u0 + (u1 - u0) * k / G : v0 + (v1 - v0) * k / G, hi = alongU ? u0 + (u1 - u0) * (k + 1) / G : v0 + (v1 - v0) * (k + 1) / G;
						if (!whole)
						{
							// on a closed surface the face may lie a turn or two away from the range
							const double least = alongU ? minU : minV, most_ = alongU ? maxU : maxV, period = alongU ? u1 - u0 : v1 - v0;
							bool covered = false;
							for (int turn = ((alongU ? s.wrapU : s.wrapV) ? -2 : 0); turn <= ((alongU ? s.wrapU : s.wrapV) ? 2 : 0) && !covered; turn++)
								covered = !(hi + turn * period < least || lo + turn * period > most_);
							if (!covered) continue;
						}
						const glm::dvec3 step = alongU ? grid[(k + 1) * (G + 1) + line] - grid[k * (G + 1) + line] : grid[line * (G + 1) + k + 1] - grid[line * (G + 1) + k];
						const double len = glm::length(step); if (len < 1e-12) continue;
						const glm::dvec3 dir = step / len;
						if (has) sum += std::acos(glm::clamp(glm::dot(dir, before), -1.0, 1.0));
						before = dir; has = true;
					}
					most = std::max(most, sum);
				}
				if (whole) return most * G / (G - 1.0);
				// the coarse picture sees a turn only between its own steps: allow for one more
				const double span = std::min(1.0, alongU ? (maxU - minU) / (u1 - u0) : (maxV - minV) / (v1 - v0));
				return most * std::min(1.0, span * G / std::max(1.0, span * G - 1.0) );
			};
			const int stepsU = glm::clamp(static_cast<int>(std::ceil(turning(true) / s.turnStep)), 1, TRIM_MAX_STEPS);
			const int stepsV = glm::clamp(static_cast<int>(std::ceil(turning(false) / s.turnStep)), 1, TRIM_MAX_STEPS);
			const double stepU = (maxU - minU) / stepsU, stepV = (maxV - minV) / stepsV;
			// Where the surface closes in (towards a pole) the points of a row lie closer together than
			// elsewhere: such a row keeps only every few of them, so that its points stay as far apart
			std::vector<int> every(stepsV + 1, 1);
			if (stepsU > 2 && stepsV > 1)
			{
				std::vector<double> apart(stepsV + 1, 0);
				double widestApart = 0;
				const double uMid = minU + stepU * (stepsU / 2);
				for (int j = 1; j < stepsV; j++)
				{
					const uv_point_t a = s.wrapped({ uMid, minV + stepV * j }), b = s.wrapped({ uMid + stepU, minV + stepV * j });
					apart[j] = glm::distance(s.point(a.x, a.y), s.point(b.x, b.y));
					widestApart = std::max(widestApart, apart[j]);
				}
				for (int j = 1; j < stepsV; j++)
					every[j] = apart[j] > 0 ? glm::clamp(static_cast<int>(widestApart / apart[j] + 0.25), 1, stepsU) : stepsU;
			}
			for (int i = 1; i < stepsU; i++)
				for (int j = 1; j < stepsV; j++)
				{
					if (i % every[j] != 0) continue;
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
					const uv_point_t in = s.wrapped({ u, v });
					positions.push_back(s.point(in.x, in.y));
				}

			// One outline, nothing inside it, and the surface straight one way: a strip between the outline's two
			// sides, each run along the direction the surface turns in
			CDT::TriangleVec triangles;
			bool strip = false;
			bool planeAboutLowPole = false;
			if (round.empty() && loops.size() == 1 && vertices.size() == boundCount && (stepsU == 1 || stepsV == 1))
			{
				for (int tryAxis = 0; tryAxis < 2 && !strip; tryAxis++)
				{
					const bool alongU = (stepsV == 1) == (tryAxis == 0);
					if (!alongU && stepsU != 1) continue;
					if (alongU && stepsV != 1) continue;
					const auto along = [&](size_t i) { return alongU ? vertices[i].x : vertices[i].y; };
					const size_t n = vertices.size();
					size_t least = 0, most = 0;
					for (size_t i = 1; i < n; i++) { if (along(i) < along(least)) least = i; if (along(i) > along(most)) most = i; }
					if (least == most) continue;
					std::vector<CDT::VertInd> one, other;
					for (size_t i = least;; i = (i + 1) % n) { one.push_back(static_cast<CDT::VertInd>(i)); if (i == most) break; }
					for (size_t i = least;; i = (i + n - 1) % n) { other.push_back(static_cast<CDT::VertInd>(i)); if (i == most) break; }
					// each side must run strictly onwards: a step straight across (a notch in the outline) would be
					// bridged by a triangle outside the face. Only the two ends of the strip may run straight across.
					bool rising = true;
					for (auto const* side : { &one, &other })
					{
						// straight across at the start, then strictly onwards, then straight across at the end
						int stage = 0;
						for (size_t i = 1; i < side->size() && rising; i++)
						{
							const double step = along((*side)[i]) - along((*side)[i - 1]);
							if (step < 0) rising = false;
							else if (step > 0) { if (stage == 2) rising = false; stage = 1; }
							else if (stage == 1) stage = 2;
						}
					}
					if (!rising) continue;
					triangles.clear();
					size_t i = 0, j = 0;
					while (i + 1 < one.size() || j + 1 < other.size())
					{
						const bool takeOne = j + 1 >= other.size() || (i + 1 < one.size() && along(one[i + 1]) <= along(other[j + 1]));
						const CDT::VertInd next = takeOne ? one[i + 1] : other[j + 1];
						if (one[i] != other[j] && next != one[i] && next != other[j])
						{
							CDT::Triangle made;
							made.vertices = { one[i], other[j], next };
							triangles.push_back(made);
						}
						if (takeOne) i++; else j++;
					}
					strip = !triangles.empty();
				}
			}
			if (!strip)
			{
			double stretchU = 1, stretchV = 1;
			// The parameters of every vertex, whatever plane the triangles are chosen in
			std::vector<uv_point_t> params(vertices.size());
			for (size_t i = 0; i < vertices.size(); i++) params[i] = { vertices[i].x, vertices[i].y };
			// A face that reaches one pole: a map centred on the pole, distance from it turned by the angle about
			// the axis. The pole is then one point, and the surface round it is laid out nearly true to size.
			const bool aboutPole = s.wrapU && reachesLow != reachesHigh;
			planeAboutLowPole = aboutPole && reachesLow;
			const glm::dvec3 tip = reachesLow ? tipLow : tipHigh;
			const auto inPlane = [&](uv_point_t const& uv, glm::dvec3 const& at)
			{
				const double far = glm::distance(at, tip);
				return CDT::V2d<double>{ far * std::cos(uv.x), far * std::sin(uv.x) };
			};
			if (aboutPole)
				for (size_t i = 0; i < vertices.size(); i++) vertices[i] = inPlane(params[i], positions[i]);
			else
			// The parameter plane stretched to the surface's own distances: there the triangulation picks the
			// triangles that are well shaped on the surface
			{
				double alongU = s.lengthU, alongV = s.lengthV;
				for (int line = 0; line <= G && !s.plain; line++)
					for (int k = 0; k < G; k++)
					{
						alongU += glm::distance(grid[(k + 1) * (G + 1) + line], grid[k * (G + 1) + line]);
						alongV += glm::distance(grid[line * (G + 1) + k + 1], grid[line * (G + 1) + k]);
					}
				// one unit is one allowed turn of the surface, or a twentieth of its longer side
				const double longer = std::max(alongU, alongV);
				const double perU = longer > 0 ? (turning(true, true) / s.turnStep + 0.05 * alongU / longer) / (u1 - u0) : 0;
				const double perV = longer > 0 ? (turning(false, true) / s.turnStep + 0.05 * alongV / longer) / (v1 - v0) : 0;
				if (perU > 0 && perV > 0 && std::isfinite(perU) && std::isfinite(perV))
				{
					for (auto& v : vertices) { v.x *= perU; v.y *= perV; }
					stretchU = perU; stretchV = perV;
				}
			}

			// Bound points that land on one another in the parameter plane (a pole, a seam met twice) are one vertex
			const auto duplicates = CDT::RemoveDuplicatesAndRemapEdges(vertices, edges);
			std::vector<char> onBound(vertices.size(), 0);
			for (size_t i = 0; i < boundCount; i++) onBound[duplicates.duplicates.empty() ? i : duplicates.mapping[i]] = 1;
			if (!duplicates.duplicates.empty())
			{
				std::vector<glm::dvec3> kept(vertices.size());
				for (size_t i = 0; i < positions.size(); i++) kept[duplicates.mapping[i]] = positions[i];
				// the first of each group wins, as the remap keeps it
				for (size_t i = positions.size(); i-- > 0;) kept[duplicates.mapping[i]] = positions[i];
				positions = std::move(kept);
				std::vector<uv_point_t> keptParams(vertices.size());
				for (size_t i = params.size(); i-- > 0;) keptParams[duplicates.mapping[i]] = params[i];
				params = std::move(keptParams);
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

			(void)boundCount;
			const auto edgeKey = [](CDT::VertInd a, CDT::VertInd b) { return (static_cast<uint64_t>(std::min(a, b)) << 32) | static_cast<uint64_t>(std::max(a, b)); };
			std::vector<uint64_t> fixed;
			fixed.reserve(once.size());
			for (auto const& e : once) fixed.push_back(edgeKey(e.v1(), e.v2()));
			std::sort(fixed.begin(), fixed.end());
			const size_t mostPoints = vertices.size() * 4 + 128;
			// where the surface itself puts each vertex (a bound point may sit a hair off it)
			std::vector<glm::dvec3> onSurface;
			std::vector<char> known;
			// The outline is triangulated once; points the refinement asks for are added to that triangulation
			CDT::Triangulation<double> whole{ CDT::VertexInsertionOrder::Auto, CDT::IntersectingConstraintEdges::NotAllowed, 0.0 };
			whole.insertVertices(vertices);
			whole.insertEdges(once);
			for (int pass = 0;; pass++)
			{
				{
					CDT::Triangulation<double> inside = whole;
					inside.eraseOuterTrianglesAndHoles();
					triangles = std::move(inside.triangles);
				}
				if (triangles.empty()) return false;
				if (pass == TRIM_REFINE_ROUNDS || vertices.size() >= mostPoints) break;
				// an inner edge that strays from the surface gets a point at its middle
				onSurface.resize(vertices.size()); known.resize(vertices.size(), 0);
				const auto surfaceAt = [&](size_t v) -> glm::dvec3 const&
				{
					if (s.even) return positions[v];
					if (!known[v])
					{
						const uv_point_t in = s.wrapped(params[v]);
						onSurface[v] = s.point(in.x, in.y); known[v] = 1;
					}
					return onSurface[v];
				};
				std::vector<CDT::V2d<double>> more;
				std::vector<glm::dvec3> moreAt;
				for (auto const& triangle : triangles)
					for (int k = 0; k < 3; k++)
					{
						const auto i = triangle.vertices[k], j = triangle.vertices[(k + 1) % 3];
						if (i > j) continue;
						if (pass >= 2 && ((i < onBound.size() && onBound[i]) || (j < onBound.size() && onBound[j]))) continue;
						if (std::binary_search(fixed.begin(), fixed.end(), edgeKey(i, j))) continue;
						const double length = glm::distance(positions[i], positions[j]);
						if (length < close * 10) continue;
						// one unit of the stretched plane is one allowed turn
						const double span2 = aboutPole ? (s.unit > 0 ? (length / s.unit) * (length / s.unit) : 1e9)
							: (vertices[i].x - vertices[j].x) * (vertices[i].x - vertices[j].x) + (vertices[i].y - vertices[j].y) * (vertices[i].y - vertices[j].y);
						if (span2 < (s.even ? TRIM_EVEN_SKIP : 0.04)) continue;
						// across a crease the split goes on the crease
						double share = 0.5;
						if (s.creased >= 0)
						{
							const double a = (s.creased == 0 ? vertices[i].x / stretchU : vertices[i].y / stretchV), b = (s.creased == 0 ? vertices[j].x / stretchU : vertices[j].y / stretchV);
							const double crease = std::round((a + b) / 2);
							if (crease > std::min(a, b) + 1e-6 && crease < std::max(a, b) - 1e-6) share = (crease - a) / (b - a);
							else
							{
								// between creases the surface is straight along that parameter: only the other counts
								const double other = s.creased == 0 ? vertices[i].y - vertices[j].y : vertices[i].x - vertices[j].x;
								if (other * other < TRIM_EVEN_SKIP) continue;
							}
						}
						// the middle in the parameters: at the pole the angle is the other end's, and the two ends are
						// taken on the same turn
						uv_point_t uvI = params[i], uvJ = params[j];
						if (aboutPole)
						{
							if (glm::distance(positions[i], tip) < close) uvI.x = uvJ.x;
							else if (glm::distance(positions[j], tip) < close) uvJ.x = uvI.x;
							uvJ.x += std::round((uvI.x - uvJ.x) / periodU) * periodU;
						}
						const uv_point_t midUV = uvI + (uvJ - uvI) * share;
						const uv_point_t in = s.wrapped(midUV);
						const glm::dvec3 on = s.point(in.x, in.y);
						const CDT::V2d<double> mid = aboutPole ? inPlane(midUV, on) : CDT::V2d<double>{ midUV.x * stretchU, midUV.y * stretchV };
						const glm::dvec3 atI = surfaceAt(i), atJ = surfaceAt(j);
						if (glm::distance(on, atI + (atJ - atI) * share) <= std::max(std::max(close * 0.5, s.allow), length * s.stray)) continue;
						more.push_back(mid); moreAt.push_back(on); params.push_back(midUV);
					}
				if (more.empty()) break;
				vertices.insert(vertices.end(), more.begin(), more.end());
				positions.insert(positions.end(), moreAt.begin(), moreAt.end());
				whole.insertVertices(more);
			}
			}
			for (auto const& triangle : triangles)
			{
				// the map about a pole at the low end of the second parameter is the parameter plane mirrored
				const bool turned = planeAboutLowPole;
				const auto i0 = triangle.vertices[0], i1 = triangle.vertices[1], i2 = triangle.vertices[2];
				const auto& a = vertices[i0]; const auto& b = vertices[i1]; const auto& c = vertices[i2];
				glm::dvec3 p0 = positions[i0], p1 = positions[i1], p2 = positions[i2];
				if (glm::length(glm::cross(p1 - p0, p2 - p0)) < 1e-18) continue;
				// Positive UV winding follows dS/du x dS/dv, regardless of bound order.
				if (((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) < 0) != turned) std::swap(p1, p2);
				geometry.AddFace(p0, p1, p2);
			}
			return true;
		}
		catch (...)
		{
			return false;
		}
	}

	bool Nurbs::fill_geometry_trimmed()
	{
		TrimSurface s;
		s.u0 = range_knots_u.x; s.u1 = range_knots_u.y; s.v0 = range_knots_v.x; s.v1 = range_knots_v.y;
		s.findWraps = true;
		s.outerTurn = this->outerTurn;
		s.turnStep = TRIM_FREEFORM_TURN_STEP;
		s.stray = TRIM_FREEFORM_STRAY;
		s.point = [this](double u, double v) { return glm::dvec3(tinynurbs::surfacePoint(*this->nurbs, u, v)); };
		s.project = [this, &s](glm::dvec3 const& pt, glm::dvec2& uv)
		{
			this->wrapU = s.wrapU; this->wrapV = s.wrapV;
			return this->project(pt, uv);
		};
		return TriangulateTrimmedSurface(this->geometry, this->bounds, s, this->scaling);
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