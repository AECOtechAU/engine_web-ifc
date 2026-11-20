/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

// Small openings in flat faces are filled, and what stands in them or behind them goes with them: bolt and
// LED holes, slots, button and lettering outlines. The rule the .frag exporter uses in Revit, with the same
// sizes: an opening under 25 mm across (or a slot under 25 mm wide and up to 1 m long).

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <unordered_map>
#include <glm/glm.hpp>
#include "../representation/geometry.h"

namespace webifc::geometry
{
	constexpr double HOLE_MAX_M = 0.025;
	constexpr double SLOT_MAX_LENGTH_M = 1.0;
	constexpr double RAISED_MAX_M = 0.100;
	constexpr double TEXT_MAX_HEIGHT_M = 0.150;

	// Off for comparing against the engine without the rule
	inline bool &FillSmallOpeningsOn()
	{
		static bool on = true;
		return on;
	}

	inline double LoopArea(const std::vector<glm::dvec2> &poly)
	{
		double area = 0;
		for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
			area += poly[j].x * poly[i].y - poly[i].x * poly[j].y;
		return area / 2;
	}

	// A hole, slot, button or letter outline: narrow (under the limit across, at most the slot length long) or
	// thin for its size (4 x area / perimeter under the limit: a circle's diameter, a letter's stroke about
	// twice). Lengths in the outline's own units; `metres` turns them into metres.
	inline bool IsSmallOpening(const std::vector<glm::dvec2> &poly, double metres)
	{
		if (poly.size() < 3) return false;
		const double max = HOLE_MAX_M / metres, slotMax = SLOT_MAX_LENGTH_M / metres;
		const double area = std::abs(LoopArea(poly));
		double perimeter = 0;
		for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) perimeter += glm::distance(poly[i], poly[j]);
		if (perimeter > 0 && 4 * area / perimeter < max) return true;
		// Narrowest width over the outline's own edge directions (a slot at any angle)
		const size_t step = std::max<size_t>(1, poly.size() / 90);
		double width = std::numeric_limits<double>::max(), length = 0;
		for (size_t i = 0; i < poly.size(); i += step)
		{
			const size_t j = (i + 1) % poly.size();
			glm::dvec2 d = poly[j] - poly[i];
			const double len = glm::length(d);
			if (len < 1e-12) continue;
			d /= len;
			double wLo = std::numeric_limits<double>::max(), wHi = -wLo, lLo = wLo, lHi = -wLo;
			for (const auto &q : poly)
			{
				const double along = q.x * d.x + q.y * d.y, across = -q.x * d.y + q.y * d.x;
				lLo = std::min(lLo, along); lHi = std::max(lHi, along);
				wLo = std::min(wLo, across); wHi = std::max(wHi, across);
			}
			if (wHi - wLo < width) { width = wHi - wLo; length = lHi - lLo; }
		}
		return width < max && length <= slotMax;
	}

	inline bool InsideOrOnLoop(const std::vector<glm::dvec2> &poly, double u, double v, double tol)
	{
		bool inside = false;
		for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
		{
			if ((poly[i].y > v) != (poly[j].y > v)
				&& u < (poly[j].x - poly[i].x) * (v - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x) inside = !inside;
		}
		if (inside) return true;
		const glm::dvec2 p(u, v);
		for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
		{
			const glm::dvec2 a = poly[j], ab = poly[i] - poly[j];
			const double l2 = glm::dot(ab, ab);
			const double t = l2 > 0 ? glm::clamp(glm::dot(p - a, ab) / l2, 0.0, 1.0) : 0.0;
			if (glm::distance(p, a + ab * t) <= tol) return true;
		}
		return false;
	}

	// A profile's small voids are filled: the solid swept from it has no small holes through it
	inline void DropSmallHoles(IfcProfile &profile, double metres)
	{
		if (!FillSmallOpeningsOn() || profile.holes.empty()) return;
		profile.holes.erase(std::remove_if(profile.holes.begin(), profile.holes.end(), [&](const IfcCurve &hole)
		{
			std::vector<glm::dvec2> poly;
			poly.reserve(hole.points.size());
			for (const auto &p : hole.points)
			{
				if (std::abs(p.z) > 1e-9) return false; // not drawn flat: left alone
				poly.emplace_back(p.x, p.y);
			}
			return IsSmallOpening(poly, metres);
		}), profile.holes.end());
	}

	// An opening cut out of an element (its mesh, in the file's units) that is small: narrow both ways across
	// (under 25 mm) and no longer than a slot (1 m). It is measured along its own faces, so it can sit at any
	// angle. Such an opening is not cut.
	inline bool IsSmallCutout(const std::vector<IfcGeometry> &pieces, double metres)
	{
		if (!FillSmallOpeningsOn()) return false;
		std::vector<glm::dvec3> points, normals;
		for (const auto &g : pieces)
		{
			if (g.halfSpace) return false;
			for (uint32_t i = 0; i < g.numPoints; i++) points.push_back(g.GetPoint(i));
			for (uint32_t i = 0; i < g.numFaces && normals.size() < 96; i++)
			{
				const auto f = g.GetFace(i);
				const glm::dvec3 n = glm::cross(g.GetPoint(f.i1) - g.GetPoint(f.i0), g.GetPoint(f.i2) - g.GetPoint(f.i0));
				const double len = glm::length(n);
				if (len < 1e-18) continue;
				const glm::dvec3 unit = n / len;
				bool known = false;
				for (const auto &m : normals)
					if (std::abs(glm::dot(m, unit)) > 0.9995) { known = true; break; }
				if (!known) normals.push_back(unit);
			}
		}
		if (points.size() < 4 || normals.empty()) return false;
		auto extent = [&](const glm::dvec3 &n)
		{
			double lo = std::numeric_limits<double>::max(), hi = -lo;
			for (const auto &p : points) { const double h = glm::dot(p, n); lo = std::min(lo, h); hi = std::max(hi, h); }
			return hi - lo;
		};
		// cheap reject on the world box first
		glm::dvec3 lo(std::numeric_limits<double>::max()), hi(-std::numeric_limits<double>::max());
		for (const auto &p : points) { lo = glm::min(lo, p); hi = glm::max(hi, p); }
		const double slot = SLOT_MAX_LENGTH_M / metres, hole = HOLE_MAX_M / metres;
		if (std::max(hi.x - lo.x, std::max(hi.y - lo.y, hi.z - lo.z)) > 1.75 * slot) return false;
		size_t first = 0;
		double e0 = std::numeric_limits<double>::max();
		for (size_t i = 0; i < normals.size(); i++) { const double e = extent(normals[i]); if (e < e0) { e0 = e; first = i; } }
		if (e0 >= hole) return false;
		glm::dvec3 second(0);
		double e1 = std::numeric_limits<double>::max();
		for (size_t i = 0; i < normals.size(); i++)
		{
			if (std::abs(glm::dot(normals[i], normals[first])) > 0.05) continue;
			const double e = extent(normals[i]);
			if (e < e1) { e1 = e; second = normals[i]; }
		}
		if (glm::length(second) < 0.5)
		{
			second = std::abs(normals[first].x) < 0.9 ? glm::dvec3(1, 0, 0) : glm::dvec3(0, 1, 0);
			second = glm::normalize(second - normals[first] * glm::dot(second, normals[first]));
			e1 = extent(second);
		}
		const double e2 = extent(glm::normalize(glm::cross(normals[first], second)));
		double e[3] = {e0, e1, e2};
		std::sort(e, e + 3);
		return e[1] < hole && e[2] <= slot;
	}

	// One face of a shell, as the small-opening rule sees it
	struct ShellFace
	{
		std::vector<IfcBound3D> bounds;
		bool planar = false;	 // flat: its openings can be filled
		bool round = false;		 // on a cylinder (axis and origin below)
		bool concave = false;	 // a cylinder seen from inside: a hole's own wall
		glm::dvec3 normal{0};	 // planar: out of the solid; zero to take it from the outline's winding
		glm::dvec3 axisOrigin{0}, axis{0};
		bool filled = false, dropped = false;
	};

	struct FilledLoop
	{
		size_t owner;
		glm::dvec3 o, x, y, n;
		std::vector<glm::dvec2> poly;
		std::vector<glm::dvec3> points;
		double minU, maxU, minV, maxV;
	};

	inline glm::dvec3 NewellNormal(const std::vector<glm::dvec3> &pts)
	{
		glm::dvec3 n(0);
		for (size_t i = 0, j = pts.size() - 1; i < pts.size(); j = i++)
		{
			n.x += (pts[j].y - pts[i].y) * (pts[j].z + pts[i].z);
			n.y += (pts[j].z - pts[i].z) * (pts[j].x + pts[i].x);
			n.z += (pts[j].x - pts[i].x) * (pts[j].y + pts[i].y);
		}
		return n;
	}

	// Fills the small openings of a shell's flat faces (their inner bounds are removed) and marks the faces that
	// sat in or behind them as dropped. Returns how many openings were filled.
	inline int FillSmallOpenings(std::vector<ShellFace> &faces, double metres)
	{
		if (!FillSmallOpeningsOn()) return 0;
		std::vector<FilledLoop> loops;
		for (size_t f = 0; f < faces.size(); f++)
		{
			auto &face = faces[f];
			if (!face.planar || face.bounds.size() < 2) continue;
			// the outline: the bound called outer, else the biggest
			size_t outer = face.bounds.size();
			for (size_t i = 0; i < face.bounds.size(); i++)
				if (face.bounds[i].type == IfcBoundType::OUTERBOUND) { outer = i; break; }
			std::vector<glm::dvec3> normals(face.bounds.size());
			if (outer == face.bounds.size())
			{
				double best = -1;
				for (size_t i = 0; i < face.bounds.size(); i++)
				{
					normals[i] = NewellNormal(face.bounds[i].curve.points);
					if (glm::length(normals[i]) > best) { best = glm::length(normals[i]); outer = i; }
				}
			}
			const auto &outline = face.bounds[outer].curve.points;
			if (outline.size() < 3) continue;
			glm::dvec3 n = NewellNormal(outline);
			if (glm::length(n) < 1e-18) continue;
			n = glm::normalize(n);
			if (glm::length(face.normal) > 0.5) n = glm::dot(n, face.normal) < 0 ? -glm::normalize(face.normal) : glm::normalize(face.normal);
			const glm::dvec3 o = outline[0];
			glm::dvec3 x = std::abs(n.x) < 0.9 ? glm::dvec3(1, 0, 0) : glm::dvec3(0, 1, 0);
			x = glm::normalize(x - n * glm::dot(x, n));
			const glm::dvec3 y = glm::cross(n, x);
			auto flat = [&](const std::vector<glm::dvec3> &pts)
			{
				std::vector<glm::dvec2> poly;
				poly.reserve(pts.size());
				for (const auto &p : pts) poly.emplace_back(glm::dot(p - o, x), glm::dot(p - o, y));
				return poly;
			};
			const auto outlinePoly = flat(outline);
			const double outlineArea = std::abs(LoopArea(outlinePoly));
			std::vector<IfcBound3D> kept;
			kept.reserve(face.bounds.size());
			for (size_t i = 0; i < face.bounds.size(); i++)
			{
				if (i == outer) { kept.push_back(face.bounds[i]); continue; }
				auto poly = flat(face.bounds[i].curve.points);
				// an opening sits inside the outline and is smaller than it
				if (poly.size() < 3 || std::abs(LoopArea(poly)) >= outlineArea || !IsSmallOpening(poly, metres))
				{
					kept.push_back(face.bounds[i]);
					continue;
				}
				FilledLoop loop;
				loop.owner = f; loop.o = o; loop.x = x; loop.y = y; loop.n = n;
				loop.minU = loop.minV = std::numeric_limits<double>::max();
				loop.maxU = loop.maxV = -std::numeric_limits<double>::max();
				for (const auto &q : poly)
				{
					loop.minU = std::min(loop.minU, q.x); loop.maxU = std::max(loop.maxU, q.x);
					loop.minV = std::min(loop.minV, q.y); loop.maxV = std::max(loop.maxV, q.y);
				}
				loop.poly = std::move(poly);
				loop.points = face.bounds[i].curve.points;
				loops.push_back(std::move(loop));
			}
			if (kept.size() != face.bounds.size())
			{
				face.bounds = std::move(kept);
				face.filled = true;
			}
		}
		if (loops.empty()) return 0;

		// Everything inside a filled opening goes with it: a hole's or a slot's walls and floor behind it, a
		// button's or a letter's sides and top standing proud of it.
		const double tol = 0.0001 / metres;
		for (size_t f = 0; f < faces.size(); f++)
		{
			auto &face = faces[f];
			if (face.filled) continue;
			glm::dvec3 c(0);
			size_t count = 0;
			glm::dvec3 lo3(std::numeric_limits<double>::max()), hi3(-std::numeric_limits<double>::max());
			for (const auto &b : face.bounds)
				for (const auto &p : b.curve.points) { c += p; count++; lo3 = glm::min(lo3, p); hi3 = glm::max(hi3, p); }
			if (count == 0) continue;
			c /= static_cast<double>(count);
			const double piece = std::max(hi3.x - lo3.x, std::max(hi3.y - lo3.y, hi3.z - lo3.z));
			for (const auto &l : loops)
			{
				if (l.owner == f) continue;
				const glm::dvec3 d = c - l.o;
				const double cu = glm::dot(d, l.x), cv = glm::dot(d, l.y);
				const double extent = std::max(l.maxU - l.minU, l.maxV - l.minV);
				// A round wall's far end is sampled on its own, so its points can sit off this outline by a
				// facet's sagitta
				const double edge = tol + 0.15 * std::min(l.maxU - l.minU, l.maxV - l.minV);
				if (cu < l.minU - edge || cu > l.maxU + edge || cv < l.minV - edge || cv > l.maxV + edge) continue;
				// Proud of the face only for small things (a button, a letter), and never further proud than
				// the opening is wide; behind it (a hole's walls) at any depth.
				const double raised = RAISED_MAX_M / metres;
				const double proud = extent <= TEXT_MAX_HEIGHT_M / metres
					? (extent < 0.008 / metres ? raised : std::min(raised, std::max(extent, 0.005 / metres)))
					: tol;
				bool all = true;
				double lo = std::numeric_limits<double>::max(), hi = -lo;
				for (const auto &b : face.bounds)
				{
					for (const auto &p : b.curve.points)
					{
						const glm::dvec3 q = p - l.o;
						const double h = glm::dot(q, l.n);
						if (h > proud || !InsideOrOnLoop(l.poly, glm::dot(q, l.x), glm::dot(q, l.y), edge)) { all = false; break; }
						lo = std::min(lo, h); hi = std::max(hi, h);
					}
					if (!all) break;
				}
				if (!all) continue;
				// Behind the opening, only what faces INTO it goes (a hole's walls and floor). A pipe running
				// on behind a bracket it passes through faces out from it and stays.
				if (lo < -tol && !face.concave && piece >= 0.25 * extent && (hi > tol || -lo > std::max(2 * extent, 0.030 / metres)))
				{
					bool out = true; // a surface this cannot judge stays
					const double mu = (l.minU + l.maxU) / 2, mv = (l.minV + l.maxV) / 2;
					if (face.planar)
					{
						glm::dvec3 fn = face.normal;
						if (glm::length(fn) < 0.5 && !face.bounds.empty()) fn = NewellNormal(face.bounds[0].curve.points);
						const double len = glm::length(fn);
						const double ru = cu - mu, rv = cv - mv, r = std::sqrt(ru * ru + rv * rv);
						out = len > 0 && r > 1e-12 && (glm::dot(fn, l.x) * ru + glm::dot(fn, l.y) * rv) / (r * len) > 0.2;
					}
					else if (face.round)
					{
						double score = 0;
						size_t n = 0;
						for (const auto &b : face.bounds)
							for (const auto &p : b.curve.points)
							{
								glm::dvec3 radial = p - face.axisOrigin;
								radial -= face.axis * glm::dot(radial, face.axis);
								const double len = glm::length(radial);
								const glm::dvec3 q = p - l.o;
								const double ru = glm::dot(q, l.x) - mu, rv = glm::dot(q, l.y) - mv, r = std::sqrt(ru * ru + rv * rv);
								if (len < 1e-12 || r < 1e-12) continue;
								score += (glm::dot(radial, l.x) * ru + glm::dot(radial, l.y) * rv) / (r * len);
								n++;
							}
						out = n > 0 && score > 0.2 * n;
					}
					if (out) continue;
				}
				face.dropped = true;
				break;
			}
		}
		// A dropped face must leave nothing open: each of its edges is shared with another dropped face or runs
		// along a filled opening. One that borders a face that stays (the opening's other end was not a hole in
		// one flat face, so it was not filled) is put back, and so on until nothing changes.
		struct EdgeUse { std::vector<uint32_t> faces; bool filledLoop = false; };
		std::unordered_map<std::string, EdgeUse> edges;
		const double grid = metres * 1e6; // points a thousandth of a millimetre apart are one point
		auto keyOf = [&](const glm::dvec3 &a, const glm::dvec3 &b)
		{
			long long q[6] = {std::llround(a.x * grid), std::llround(a.y * grid), std::llround(a.z * grid),
							  std::llround(b.x * grid), std::llround(b.y * grid), std::llround(b.z * grid)};
			if (std::lexicographical_compare(q + 3, q + 6, q, q + 3)) { std::swap(q[0], q[3]); std::swap(q[1], q[4]); std::swap(q[2], q[5]); }
			return std::string(reinterpret_cast<const char *>(q), sizeof(q));
		};
		for (size_t f = 0; f < faces.size(); f++)
			for (const auto &bound : faces[f].bounds)
			{
				const auto &pts = bound.curve.points;
				for (size_t i = 0, j = pts.size() - 1; i < pts.size(); j = i++)
					if (pts[i] != pts[j]) edges[keyOf(pts[j], pts[i])].faces.push_back(static_cast<uint32_t>(f));
			}
		for (const auto &l : loops)
			for (size_t i = 0, j = l.points.size() - 1; i < l.points.size(); j = i++)
			{
				auto it = edges.find(keyOf(l.points[j], l.points[i]));
				if (it != edges.end()) it->second.filledLoop = true;
			}
		for (bool changed = true; changed;)
		{
			changed = false;
			for (size_t f = 0; f < faces.size(); f++)
			{
				if (!faces[f].dropped) continue;
				bool open = false;
				for (const auto &bound : faces[f].bounds)
				{
					const auto &pts = bound.curve.points;
					for (size_t i = 0, j = pts.size() - 1; i < pts.size() && !open; j = i++)
					{
						if (pts[i] == pts[j]) continue;
						const auto &use = edges[keyOf(pts[j], pts[i])];
						if (use.filledLoop) continue;
						for (const auto other : use.faces)
							if (other != f && !faces[other].dropped) { open = true; break; }
					}
					if (open) break;
				}
				if (open) { faces[f].dropped = false; changed = true; }
			}
		}
		return static_cast<int>(loops.size());
	}
}
