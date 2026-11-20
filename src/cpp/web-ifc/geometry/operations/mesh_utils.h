/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

#pragma once

#include "../nurbs.h"
#include <vector>
#include <array>
#include <unordered_map>
#include <optional>
#include <memory>
#include <limits>
#include <cmath>
#include <cstdint>
#include <tinynurbs/tinynurbs.h>
#include <spdlog/spdlog.h>
#include "geometryutils.h"
#include <glm/glm.hpp>

#define CONST_PI 3.141592653589793238462643383279502884L

namespace webifc::geometry
{

	// TODO: This function is wrong, but changing it breaks many models #172 -> 15.ifc, should be reviewed
	inline double VectorToAngle3D(double x, double y)
	{
		double dd = sqrt(x * x + y * y);
		if (std::abs(dd) < EPS_MINISCULE) {
			return 0;
		}
		double xx = x / dd;
		double yy = y / dd;

		double angle = asin(xx);
		double cosv = cos(angle);

		if (glm::abs(yy - cosv) > 1e-5)
		{
			angle = acos(yy);
			double sinv = sin(angle);
			cosv = cos(angle);
			if (glm::abs(yy - cosv) > 1e-5 || glm::abs(xx - sinv) > 1e-5)
			{
				angle = angle + (CONST_PI - angle) * 2;
				sinv = sin(angle);
				cosv = cos(angle);
				if (glm::abs(yy - cosv) > 1e-5 || glm::abs(xx - sinv) > 1e-5)
				{
					angle = angle + CONST_PI;
				}
			}
		}

		return (angle / (2 * CONST_PI)) * 360;
	}

	// TODO: review and simplify
	inline void TriangulateRevolution(IfcGeometry& geometry, std::vector<IfcBound3D> const& bounds, IfcSurface const& surface, double numRots)
	{
		spdlog::debug("[TriangulateRevolution({})]");
			// First we get the revolution data

		glm::dmat4 transform;
		transform[3] = surface.RevolutionSurface.Direction[3];
		transform[0] = glm::normalize(surface.RevolutionSurface.Direction[0]);
		transform[1] = glm::normalize(surface.RevolutionSurface.Direction[1]);
		transform[2] = glm::normalize(surface.RevolutionSurface.Direction[2]);

		// Find the max. curve index in the boundary
		int maxTeam = 0;
		for (size_t i = 0; i < bounds.size(); i++)
		{
			for (size_t j = 0; j < bounds[i].curve.indices.size(); j++)
			{
				if (bounds[i].curve.indices[j] > maxTeam)
				{
					maxTeam = bounds[i].curve.indices[j];
				}
			}
		}

		// We group each point with their boundary
		std::vector<std::vector<glm::dvec3>> boundingGroups;

		for (int r = 0; r <= maxTeam; r++)
		{
			std::vector<glm::dvec3> boundingTemp = std::vector<glm::dvec3>();
			for (size_t i = 0; i < bounds.size(); i++)
			{
				for (size_t j = 0; j < bounds[i].curve.points.size(); j++)
				{
					if (bounds[i].curve.indices[j] == r)
					{
						boundingTemp.push_back(bounds[i].curve.points[j]);
					}
				}
			}
			boundingGroups.push_back(boundingTemp);
		}

		std::vector<glm::dvec3> bounding;
		std::vector<double> angleVec;
		std::vector<double> angleDsp;

		// Now we construct the bounding box of the boundary ...
		// ... by adding the middle point of all curves

		for (size_t i = 0; i < boundingGroups.size(); i++)
		{
			double xx = 0;
			double yy = 0;
			double zz = 0;
			double cc = 0;
			for (size_t j = 0; j < boundingGroups[i].size(); j++)
			{
				// If it is the first point of the group we close the previous group ...
				//  ... and create a new one. Else, the point is of the current group
				if (j == boundingGroups[i].size() - 1)
				{
					if (cc > 0)
					{
						xx /= cc;
						yy /= cc;
						zz /= cc;
						bounding.push_back(glm::dvec3(xx, yy, zz));
					}
					xx = boundingGroups[i][j].x;
					yy = boundingGroups[i][j].y;
					zz = boundingGroups[i][j].z;
					cc = 1;
				}
				else
				{
					xx += boundingGroups[i][j].x;
					yy += boundingGroups[i][j].y;
					zz += boundingGroups[i][j].z;
					cc++;
				}
			}
		}

		// There is a problem when points in the revolution are around 0 degrees
		// Numerical instabilities can make these points to jump from 0 to 360
		// It causes lots of trouble when drawing the boundaries in the revolution

		// The method presented here finds the angle of each point, measures the ...
		//  ... angular difference and then, if the difference is bigger than 180 ...
		//  ... corrects it to a lesser value. Finally it gets the first angle and ...
		//  ... adds the angular differences again, reconstructing a corrected boundary.

		// Now we find the angle of each point in the reference plane of the cylinder

		for (size_t j = 0; j < bounding.size(); j++)
		{
			double xx = bounding[j].x - transform[3].x;
			double yy = bounding[j].y - transform[3].y;
			double zz = bounding[j].z - transform[3].z;
			double dx = transform[0].x * xx + transform[0].y * yy + transform[0].z * zz;
			double dy = transform[1].x * xx + transform[1].y * yy + transform[1].z * zz;
				//				double dz = vecZ.x * xx + vecZ.y * yy + vecZ.z * zz;
			double temp = VectorToAngle3D(dx, dy);
			while (temp < 0)
			{
				temp += 360;
			}
			while (temp > 360)
			{
				temp -= 360;
			}
			angleVec.push_back(temp);
		}

		for (size_t i = 0; i < angleVec.size() - 1; i++)
		{
			if (angleVec[i] - angleVec[i + 1] > 180)
			{
				angleDsp.push_back(360 - (angleVec[i] - angleVec[i + 1]));
			}
			else if (angleVec[i] - angleVec[i + 1] < -180)
			{
				angleDsp.push_back(-(angleVec[i] - angleVec[i + 1] + 360));
			}
			else
			{
				angleDsp.push_back(angleVec[i + 1] - angleVec[i]);
			}
		}

		double startDegrees = angleVec[0];
		double endDegrees = angleVec[0];

		// Add angular differences starting from the first angle. We also correct the start and end angles

		double temp = angleVec[0];
		for (size_t i = 0; i < angleDsp.size(); i++)
		{
			temp += angleDsp[i];
			if (endDegrees < temp)
			{
				endDegrees = temp;
			}
			if (startDegrees > temp)
			{
				startDegrees = temp;
			}
		}

		// Canviar i passar la funció de calcular els angles inicial i final aqui de manera que sigui més simple.

		bimGeometry::Geometry geom = bimGeometry::Revolution(
			transform, 
			startDegrees,
			endDegrees,
			surface.RevolutionSurface.Profile.curve.points,
			numRots
			);
		
		for (int r = 0; r < geom.numFaces; r++)
		{
		 	bimGeometry::Face f = geom.GetFace(r);
			geometry.AddFace(geom.GetPoint(f.i0), geom.GetPoint(f.i1),  geom.GetPoint(f.i2));
		}
	}

	// TODO: review and simplify
	inline void TriangulateCylindricalSurface(IfcGeometry &geometry, std::vector<IfcBound3D> const& bounds, IfcSurface const& surface, double numRots)
	{
		spdlog::debug("[TriangulateCylindricalSurface({})]");

		glm::dvec3 cent = surface.transformation[3];
		glm::dvec3 vecX = glm::normalize(surface.transformation[0]);
		glm::dvec3 vecY = glm::normalize(surface.transformation[1]);
		glm::dvec3 vecZ = glm::normalize(surface.transformation[2]);

		// Find the max. curve index in the boundary
		int maxTeam = 0;
		for (size_t i = 0; i < bounds.size(); i++)
		{
			for (size_t j = 0; j < bounds[i].curve.indices.size(); j++)
			{
				if (bounds[i].curve.indices[j] > maxTeam)
				{
					maxTeam = bounds[i].curve.indices[j];
				}
			}
		}

		// We group each point with their boundary
		std::vector<std::vector<glm::dvec3>> boundingGroups;

		for (int r = 0; r <= maxTeam; r++)
		{
			std::vector<glm::dvec3> boundingTemp = std::vector<glm::dvec3>();
			for (size_t i = 0; i < bounds.size(); i++)
			{
				for (size_t j = 0; j < bounds[i].curve.points.size(); j++)
				{
					if (bounds[i].curve.indices[j] == r)
					{
						boundingTemp.push_back(bounds[i].curve.points[j]);
					}
				}
			}
			boundingGroups.push_back(boundingTemp);
		}

		std::vector<double> angleVec;
		std::vector<double> angleDsp;

		double minZ = 1e+10;
		double maxZ = -1e+10;

		// Find the relative coordinates of each curve point in the cylinder reference plane
		// Only retain the max and min relative Z

		for (size_t i = 0; i < boundingGroups.size(); i++)
		{
			for (size_t j = 0; j < boundingGroups[i].size(); j++)
			{
				glm::dvec3 vv = boundingGroups[i][j] - cent;
				double dz = glm::dot(vecZ, vv);
				if (maxZ < dz)
				{
					maxZ = dz;
				}
				if (minZ > dz)
				{
					minZ = dz;
				}
			}
		}

		// First we get the cylinder data
		std::vector<glm::dvec3> bounding;
		int repeats = 0;
		bool start = false;
		size_t id = 0;

		// In the case of boundary lines having only 2 endings...
		//... we omit these lines and add solely curves having > 2 points...
		//... starting from a 2 point line, by doing it this way we don't have repeated points

		while (repeats < boundingGroups.size() * 3)
		{
			if (id >= boundingGroups.size())
			{
				id = 0;
			}
			if (boundingGroups[id].size() < 3)
			{
				if (!start)
				{
					start = true;
				}
				else
				{
					break;
				}
			}
			if (boundingGroups[id].size() > 2 && start)
			{
				for (size_t i = 0; i < boundingGroups[id].size(); i++)
				{
					bounding.push_back(boundingGroups[id][i]);
				}
			}
			id++;
			repeats++;
		}

		// If the previous method finds nothing, then we don't have straight lines ...
		//... then we add all boundary points directly
		if (bounding.size() == 0)
		{
			for (size_t j = 0; j < boundingGroups.size(); j++)
			{
				for (size_t i = 0; i < boundingGroups[j].size(); i++)
				{
					bounding.push_back(boundingGroups[j][i]);
				}
			}
		}

		double startDegrees = 0;
		double endDegrees = 360;

		// Now we project the points in the cylinder surface
		// There is a problem when points in the cylinder are around 0 degrees
		// Numerical instabilities can make these points to jump from 0 to 360
		// It causes lots of trouble when drawing the boundaries in the cylinder

		// The method presented here finds the angle of each point, measures the ...
		//  ... angular difference and then, if the difference is bigger than 180 ...
		//  ... corrects it to a lesser value. Finally it gets the first angle and ...
		//  ... adds the angular differences again, reconstructing a corrected boundary.

		// Now we find the angle of each point in the reference plane of the cylinder

		for (size_t j = 0; j < bounding.size(); j++)
		{
			glm::dvec3 vv = bounding[j] - cent;
			double dx = glm::dot(vecX, vv);
			double dy = glm::dot(vecY, vv);
				// double dz = glm::dot(vecZ, vv);
			double temp = VectorToAngle3D(dx, dy);
			while (temp < 0)
			{
				temp += 360;
			}
			while (temp > 360)
			{
				temp -= 360;
			}
			angleVec.push_back(temp);
		}

		// Then we find the angular difference
		if (angleVec.size() == 0) return;

		for (size_t i = 0; i < angleVec.size() - 1; i++)
		{
			if (angleVec[i] - angleVec[i + 1] > 180)
			{
				angleDsp.push_back(360 - (angleVec[i] - angleVec[i + 1]));
			}
			else if (angleVec[i] - angleVec[i + 1] < -180)
			{
				angleDsp.push_back(-(angleVec[i] - angleVec[i + 1] + 360));
			}
			else
			{
				angleDsp.push_back(angleVec[i + 1] - angleVec[i]);
			}
		}

		startDegrees = angleVec[0];
		endDegrees = angleVec[0];

		// Add angular differences starting from the first angle. We also correct the start and end angles

		double temp = angleVec[0];
		for (size_t i = 0; i < angleDsp.size(); i++)
		{
			temp += angleDsp[i];
			if (endDegrees < temp)
			{
				endDegrees = temp;
			}
			if (startDegrees > temp)
			{
				startDegrees = temp;
			}
		}

			// Then we use the start and end angles as bounding boxes of the boundary ...
			//  ... we will represent this bounding box.

		while (startDegrees < -360)
		{
			startDegrees += 360;
		}

		bimGeometry::Geometry geom = bimGeometry::RevolveCylinder(
			surface.transformation, 
			startDegrees,
			endDegrees,
			minZ,
			maxZ,
			numRots, 
			surface.CylinderSurface.Radius);
		
		for (int r = 0; r < geom.numFaces; r++)
		{
		 	bimGeometry::Face f = geom.GetFace(r);
			geometry.AddFace(geom.GetPoint(f.i0), geom.GetPoint(f.i1),  geom.GetPoint(f.i2));
		}
	}

    // TODO: review and simplify
    inline void TriangulateExtrusion(IfcGeometry &geometry, std::vector<IfcBound3D> const& bounds, IfcSurface const& surface)
	{
		spdlog::debug("[TriangulateExtrusion({})]");

		double len = surface.ExtrusionSurface.Length;
		glm::dvec3 dir = surface.ExtrusionSurface.Direction;

		if (!surface.ExtrusionSurface.Profile.isComposite)
		{
			std::vector<glm::dvec3> points = surface.ExtrusionSurface.Profile.curve.points;
			bimGeometry::Geometry geom = bimGeometry::Extrude(points, dir, len);
			for (int r = 0; r < geom.numFaces; r++)
			{
				bimGeometry::Face f = geom.GetFace(r);
				geometry.AddFace(geom.GetPoint(f.i0), geom.GetPoint(f.i1),  geom.GetPoint(f.i2));
			}
		}
		else
		{
			for (size_t i = 0; i < surface.ExtrusionSurface.Profile.profiles.size(); i++)
			{
				std::vector<glm::dvec3> points = surface.ExtrusionSurface.Profile.profiles[i].curve.points;
                bimGeometry::Geometry geom = bimGeometry::Extrude(points, dir, len);
				for (int r = 0; r < geom.numFaces; r++)
				{
					bimGeometry::Face f = geom.GetFace(r);
					geometry.AddFace(geom.GetPoint(f.i0), geom.GetPoint(f.i1),  geom.GetPoint(f.i2));
				}
            }
		}
    }

	// A run of points taken as a curve whose parameter counts the points: 0 at the first, 1 at the second
	struct TrimPolyline
	{
		std::vector<glm::dvec2> flat;
		// The parameter of the curve point nearest a position in the plane, and how far off it is
		double nearest(glm::dvec2 const& q, double& distance) const
		{
			double best = std::numeric_limits<double>::max(), at = 0;
			for (size_t i = 0; i + 1 < flat.size(); i++)
			{
				const glm::dvec2 ab = flat[i + 1] - flat[i];
				const double len2 = glm::dot(ab, ab);
				const double t = len2 > 0 ? glm::clamp(glm::dot(q - flat[i], ab) / len2, 0.0, 1.0) : 0.0;
				const double d = glm::distance(q, flat[i] + ab * t);
				if (d < best) { best = d; at = static_cast<double>(i) + t; }
			}
			distance = best;
			return at;
		}
		template <typename T> static T at(std::vector<T> const& values, double t)
		{
			const double last = static_cast<double>(values.size() - 1);
			t = glm::clamp(t, 0.0, last);
			const size_t i = std::min(static_cast<size_t>(t), values.size() - 2);
			const double f = t - static_cast<double>(i);
			return values[i] * (1.0 - f) + values[i + 1] * f;
		}
	};

	// A smooth curve through a run of points, its parameter counting the points as TrimPolyline's does. Where the
	// run turns by more than 35 degrees at a point the curve has a corner there.
	template <typename T> struct TrimCurve
	{
		std::vector<T> points, arriving, leaving;
		void build(std::vector<T> const& run, bool closed)
		{
			points = run;
			const size_t n = points.size();
			arriving.assign(n, T(0)); leaving.assign(n, T(0));
			const auto unit = [](T const& v) { const double l = glm::length(v); return l > 0 ? v / l : v; };
			const size_t last = n - 1;
			for (size_t i = 0; i < n; i++)
			{
				const bool hasBefore = i > 0 || closed, hasAfter = i < last || closed;
				const T before = hasBefore ? points[i] - points[i > 0 ? i - 1 : last - 1] : T(0);
				const T after = hasAfter ? points[i < last ? i + 1 : 1] - points[i] : T(0);
				const double lb = glm::length(before), la = glm::length(after);
				if (lb > 0 && la > 0)
				{
					const T b = before / lb, a = after / la;
					if (glm::dot(a, b) < 0.819) { arriving[i] = b; leaving[i] = a; }
					else arriving[i] = leaving[i] = unit(b * la + a * lb);
				}
				else arriving[i] = leaving[i] = lb > 0 ? before / lb : (la > 0 ? after / la : T(0));
			}
			if (!closed && n >= 3)
			{
				// at an open end the curve keeps the bend it has at the next point
				leaving[0] = arriving[0] = unit(unit(points[1] - points[0]) * 2.0 - arriving[1]);
				arriving[last] = leaving[last] = unit(unit(points[last] - points[last - 1]) * 2.0 - leaving[last - 1]);
			}
		}
		T at(double t) const
		{
			const double end = static_cast<double>(points.size() - 1);
			t = glm::clamp(t, 0.0, end);
			const size_t i = std::min(static_cast<size_t>(t), points.size() - 2);
			const double f = t - static_cast<double>(i), f2 = f * f, f3 = f2 * f;
			const double l = glm::length(points[i + 1] - points[i]);
			return points[i] * (2 * f3 - 3 * f2 + 1) + leaving[i] * (l * (f3 - 2 * f2 + f)) + points[i + 1] * (3 * f2 - 2 * f3) + arriving[i + 1] * (l * (f3 - f2));
		}
	};

	// The parameter of the point of a smooth curve nearest a position, searched a point either side of a first
	// guess (the nearest on the straight pieces), and how far off it is
	inline double TrimNearestOnCurve(TrimCurve<glm::dvec2> const& curve, glm::dvec2 const& q, double guess, double& distance)
	{
		double lo = std::max(0.0, guess - 1), hi = std::min(static_cast<double>(curve.points.size() - 1), guess + 1);
		for (int i = 0; i < 14; i++)
		{
			const double a = lo + (hi - lo) / 3, b = hi - (hi - lo) / 3;
			if (glm::distance(curve.at(a), q) < glm::distance(curve.at(b), q)) hi = b; else lo = a;
		}
		const double t = (lo + hi) / 2;
		distance = glm::distance(curve.at(t), q);
		return t;
	}

	// S(u, v) = centre + radius * (cos u * X + sin u * Y) + v * Z
	inline bool TriangulateTrimmedCylinder(IfcGeometry &geometry, std::vector<IfcBound3D> const& bounds, IfcSurface const& surface, double scaling, int outerTurn = 0)
	{
		const glm::dvec3 centre = surface.transformation[3];
		const glm::dvec3 X = glm::normalize(glm::dvec3(surface.transformation[0])), Y = glm::normalize(glm::dvec3(surface.transformation[1])), Z = glm::normalize(glm::dvec3(surface.transformation[2]));
		const double radius = surface.CylinderSurface.Radius;
		if (!(radius > 0)) return false;
		double low = std::numeric_limits<double>::max(), high = std::numeric_limits<double>::lowest();
		for (auto const& bound : bounds)
			for (auto const& pt : bound.curve.points) { const double h = glm::dot(pt - centre, Z); low = std::min(low, h); high = std::max(high, h); }
		if (!(high > low)) return false;
		TrimSurface s;
		s.u0 = 0; s.u1 = 2 * CONST_PI; s.v0 = low; s.v1 = high; s.wrapU = true;
		s.point = [=](double u, double v) { return centre + radius * (std::cos(u) * X + std::sin(u) * Y) + v * Z; };
		s.project = [=](glm::dvec3 const& pt, glm::dvec2& uv)
		{
			const glm::dvec3 d = pt - centre;
			const double x = glm::dot(d, X), y = glm::dot(d, Y), r = std::hypot(x, y);
			double angle = r > radius * 1e-9 ? std::atan2(y, x) : uv.x;
			angle += std::round((uv.x - angle) / (2 * CONST_PI)) * (2 * CONST_PI);
			uv = { angle, glm::clamp(glm::dot(d, Z), low, high) };
			return std::abs(r - radius);
		};
		s.direct = true;
		s.outerTurn = outerTurn;
		s.even = true;
		s.plain = true;
		// meshed to the step its outline is drawn with
		s.turnStep = RoundStep(radius * scaling);
		s.stray = 1.5 * (1 - std::cos(s.turnStep / 2)) / (2 * std::sin(s.turnStep / 2));
		s.turnU = 2 * CONST_PI; s.lengthU = 2 * CONST_PI * radius; s.lengthV = high - low;
		return TriangulateTrimmedSurface(geometry, bounds, s, scaling);
	}

	// Whether a surface is a sphere, and which: one given as such, or an arc of a circle centred on the axis
	// it is turned about
	// looks: +1 when the surface's own normal points out of the ball, -1 when into it (that depends on which way
	// the profile of a turned surface runs)
	inline bool SphereOf(IfcSurface const& surface, glm::dvec3& centre, double& radius, int* looks = nullptr)
	{
		if (looks) *looks = 1;
		if (surface.SphereSurface.Active)
		{
			centre = surface.transformation[3];
			radius = surface.SphereSurface.Radius;
			return radius > 0;
		}
		if (!surface.RevolutionSurface.Active || surface.RevolutionSurface.Profile.isComposite) return false;
		auto const& given = surface.RevolutionSurface.Profile.curve.points;
		if (given.size() < 3) return false;
		const glm::dmat4 axis = surface.transformation * surface.RevolutionSurface.Direction;
		const glm::dvec3 origin = axis[3];
		const glm::dvec3 Z = glm::normalize(glm::dvec3(axis[2]));
		std::vector<glm::dvec2> flat;
		for (auto const& pt : given)
		{
			const glm::dvec3 d = glm::dvec3(surface.transformation * glm::dvec4(pt, 1.0)) - origin;
			const double h = glm::dot(d, Z);
			flat.emplace_back(glm::length(d - Z * h), h);
		}
		const glm::dvec2 a = flat.front(), b = flat[flat.size() / 2], c = flat.back();
		const double det = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
		const double reach = glm::distance(a, c);
		if (std::abs(det) < 1e-12 * reach * reach || reach <= 0) return false;
		const glm::dvec2 middle((glm::dot(a, a) * (b.y - c.y) + glm::dot(b, b) * (c.y - a.y) + glm::dot(c, c) * (a.y - b.y)) / det,
		                        (glm::dot(a, a) * (c.x - b.x) + glm::dot(b, b) * (a.x - c.x) + glm::dot(c, c) * (b.x - a.x)) / det);
		radius = glm::distance(a, middle);
		if (std::abs(middle.x) > radius * 1e-3) return false;
		for (auto const& pt : flat) if (std::abs(glm::distance(pt, middle) - radius) > radius * 2e-3) return false;
		centre = origin + Z * middle.y;
		if (looks)
		{
			// the normal is (the way round the axis) x (the way along the profile)
			const size_t i = flat.size() / 2;
			const glm::dvec2 step = flat[i] - flat[i - 1], at = (flat[i] + flat[i - 1]) * 0.5;
			*looks = step.y * at.x - step.x * (at.y - middle.y) >= 0 ? 1 : -1;
		}
		return radius > 0;
	}

	// S(u, v) = the profile point at v turned by u about the axis
	inline bool TriangulateTrimmedRevolution(IfcGeometry &geometry, std::vector<IfcBound3D> const& bounds, IfcSurface const& surface, double scaling, int outerTurn = 0)
	{
		auto const& given = surface.RevolutionSurface.Profile.curve.points;
		if (given.size() < 2 || surface.RevolutionSurface.Profile.isComposite) return false;
		// The profile and the axis are given in the frame of the surface's position
		const glm::dmat4 axis = surface.transformation * surface.RevolutionSurface.Direction;
		const glm::dvec3 origin = axis[3];
		const glm::dvec3 Z = glm::normalize(glm::dvec3(axis[2]));
		std::vector<glm::dvec3> profile(given.size());
		for (size_t i = 0; i < given.size(); i++) profile[i] = glm::dvec3(surface.transformation * glm::dvec4(given[i], 1.0));
		// each profile point by its distance from the axis, its height along it and its angle about it
		glm::dvec3 X(0); double widest = 0;
		for (auto const& pt : profile)
		{
			const glm::dvec3 d = pt - origin, across = d - Z * glm::dot(d, Z);
			if (glm::length(across) > widest) { widest = glm::length(across); X = across / widest; }
		}
		if (!(widest > 0)) return false;
		const glm::dvec3 Y = glm::cross(Z, X);
		auto line = std::make_shared<TrimPolyline>();
		auto angles = std::make_shared<std::vector<double>>();
		for (auto const& pt : profile)
		{
			const glm::dvec3 d = pt - origin;
			const double x = glm::dot(d, X), y = glm::dot(d, Y), r = std::hypot(x, y);
			line->flat.emplace_back(r, glm::dot(d, Z));
			angles->push_back(r > widest * 1e-9 ? std::atan2(y, x) : 0.0);
		}
		// a profile point on the axis has no angle of its own: it takes its neighbour's
		for (size_t i = 0; i < profile.size(); i++)
			if (line->flat[i].x <= widest * 1e-9) (*angles)[i] = (*angles)[i > 0 ? i - 1 : 1];
		// A profile that is a straight line or an arc of a circle is carried on past its ends as far as the
		// face's bounds reach, never across the axis: files give only part of the profile of some faces
		{
			auto& flat = line->flat;
			const size_t n = flat.size();
			bool oneAngle = true;
			for (double a : *angles) if (std::abs(a - angles->front()) > 1e-6) oneAngle = false;
			const double reach = glm::distance(flat.front(), flat.back());
			if (oneAngle && reach > widest * 1e-9)
			{
				const double phi = angles->front();
				std::vector<glm::dvec2> seen;
				for (auto const& bound : bounds)
					for (auto const& pt : bound.curve.points)
					{
						const glm::dvec3 d = pt - origin;
						seen.emplace_back(std::hypot(glm::dot(d, X), glm::dot(d, Y)), glm::dot(d, Z));
					}
				const glm::dvec2 along = (flat.back() - flat.front()) / reach;
				double bend = 0;
				for (auto const& pt : flat) { const glm::dvec2 d = pt - flat.front(); bend = std::max(bend, std::abs(d.x * along.y - d.y * along.x)); }
				std::vector<glm::dvec2> before, after;
				if (bend < reach * 1e-6)
				{
					// a straight line
					double lo = 0, hi = reach;
					for (auto const& q : seen)
					{
						const glm::dvec2 d = q - flat.front();
						if (std::abs(d.x * along.y - d.y * along.x) > reach * 1e-3 + widest * 1e-4) continue;
						lo = std::min(lo, glm::dot(d, along)); hi = std::max(hi, glm::dot(d, along));
					}
					// not across the axis
					if (along.x > 0) lo = std::max(lo, -flat.front().x / along.x);
					if (along.x < 0) hi = std::min(hi, -flat.front().x / along.x);
					if (lo < -reach * 1e-6) before.push_back(flat.front() + along * lo);
					if (hi > reach * (1 + 1e-6)) after.push_back(flat.front() + along * hi);
				}
				else if (n >= 3)
				{
					// an arc of a circle: its centre from the first, middle and last points
					const glm::dvec2 a = flat.front(), b = flat[n / 2], c = flat.back();
					const double det = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
					if (std::abs(det) > 1e-12 * reach * reach)
					{
						const glm::dvec2 centre((glm::dot(a, a) * (b.y - c.y) + glm::dot(b, b) * (c.y - a.y) + glm::dot(c, c) * (a.y - b.y)) / det,
						                        (glm::dot(a, a) * (c.x - b.x) + glm::dot(b, b) * (a.x - c.x) + glm::dot(c, c) * (b.x - a.x)) / det);
						const double radius = glm::distance(a, centre);
						bool round = true;
						for (auto const& pt : flat) if (std::abs(glm::distance(pt, centre) - radius) > radius * 2e-3) round = false;
						const auto turnOf = [&](glm::dvec2 const& q) { return std::atan2(q.y - centre.y, q.x - centre.x); };
						const auto toPi = [](double x) { while (x > CONST_PI) x -= 2 * CONST_PI; while (x < -CONST_PI) x += 2 * CONST_PI; return x; };
						if (round)
						{
							const double start = turnOf(a);
							const double way = toPi(turnOf(flat[1]) - start) >= 0 ? 1.0 : -1.0;
							double sweep = 0;
							for (size_t i = 1; i < n; i++) sweep += way * toPi(turnOf(flat[i]) - turnOf(flat[i - 1]));
							const double step = sweep / static_cast<double>(n - 1);
							if (sweep > 0 && step > 0)
							{
								double lo = 0, hi = sweep;
								for (auto const& q : seen)
								{
									if (std::abs(glm::distance(q, centre) - radius) > radius * 0.02) continue;
									const double at = sweep / 2 + way * toPi(turnOf(q) - (start + way * sweep / 2));
									lo = std::min(lo, at); hi = std::max(hi, at);
								}
								const auto on = [&](double at) { return centre + radius * glm::dvec2(std::cos(start + way * at), std::sin(start + way * at)); };
								// outwards from an end by the profile's own step, stopping on the axis
								const auto carry = [&](double from, double to, std::vector<glm::dvec2>& out)
								{
									const double dir = to > from ? 1.0 : -1.0;
									double at = from;
									for (int guard = 0; guard < 4096 && dir * (to - at) > step * 0.25; guard++)
									{
										double next = at + dir * step;
										if (dir * (next - to) > 0) next = to;
										if (on(next).x < 0)
										{
											double in = at, outside = next;
											for (int i = 0; i < 40; i++) { const double mid = (in + outside) / 2; if (on(mid).x < 0) outside = mid; else in = mid; }
											glm::dvec2 tip = on(in); tip.x = 0;
											out.push_back(tip);
											return;
										}
										out.push_back(on(next));
										at = next;
									}
								};
								lo = std::max(lo, sweep - static_cast<double>(2 * CONST_PI) + step);
								hi = std::min(hi, lo + static_cast<double>(2 * CONST_PI) - step);
								if (flat.front().x > widest * 1e-9) carry(0, lo, before);
								if (flat.back().x > widest * 1e-9) carry(sweep, hi, after);
							}
						}
					}
				}
				if (!before.empty() || !after.empty())
				{
					std::reverse(before.begin(), before.end());
					flat.insert(flat.begin(), before.begin(), before.end());
					flat.insert(flat.end(), after.begin(), after.end());
					angles->assign(flat.size(), phi);
					for (auto const& pt : flat) widest = std::max(widest, pt.x);
				}
			}
		}
		const double last = static_cast<double>(line->flat.size() - 1);
		TrimSurface s;
		s.u0 = 0; s.u1 = 2 * CONST_PI; s.v0 = 0; s.v1 = last; s.wrapU = true;
		const bool closed = line->flat.size() > 2 && glm::distance(line->flat.front(), line->flat.back()) < widest * 1e-9;
		s.wrapV = closed;
		// meshed to the step its outline is drawn with
		s.turnStep = RoundStep(widest * scaling);
		s.stray = 1.5 * (1 - std::cos(s.turnStep / 2)) / (2 * std::sin(s.turnStep / 2));
		// The allowance of a round thing as tight as the surface gets: its profile's tightest bend (corners
		// apart), or its widest circle about the axis if that is tighter
		{
			double tightest = widest;
			auto const& flat = line->flat;
			for (size_t i = 1; i + 1 < flat.size(); i++)
			{
				const glm::dvec2 a = flat[i] - flat[i - 1], b = flat[i + 1] - flat[i];
				const double la = glm::length(a), lb = glm::length(b);
				if (la <= 0 || lb <= 0) continue;
				const double turn = std::acos(glm::clamp(glm::dot(a, b) / (la * lb), -1.0, 1.0));
				if (turn > 1e-3 && turn < 0.61) tightest = std::min(tightest, (la + lb) / 2 / turn);
			}
			s.allow = tightest * (1 - std::cos(RoundStep(tightest * scaling) / 2));
			s.unit = 2 * tightest * std::sin(RoundStep(tightest * scaling) / 2);
		}
		auto curve = std::make_shared<TrimCurve<glm::dvec2>>();
		curve->build(line->flat, closed);
		s.point = [=](double u, double v)
		{
			const glm::dvec2 rh = curve->at(v);
			const double angle = u + TrimPolyline::at(*angles, v);
			return origin + rh.x * (std::cos(angle) * X + std::sin(angle) * Y) + rh.y * Z;
		};
		s.project = [=](glm::dvec3 const& pt, glm::dvec2& uv)
		{
			const glm::dvec3 d = pt - origin;
			const double x = glm::dot(d, X), y = glm::dot(d, Y), r = std::hypot(x, y);
			double left = 0;
			double v = line->nearest({ r, glm::dot(d, Z) }, left);
			v = TrimNearestOnCurve(*curve, { r, glm::dot(d, Z) }, v, left);
			double angle = r > widest * 1e-9 ? std::atan2(y, x) - TrimPolyline::at(*angles, v) : uv.x;
			angle += std::round((uv.x - angle) / (2 * CONST_PI)) * (2 * CONST_PI);
			if (closed) v += std::round((uv.y - v) / last) * last;
			uv = { angle, v };
			return left;
		};
		s.direct = true;
		s.outerTurn = outerTurn;
		return TriangulateTrimmedSurface(geometry, bounds, s, scaling);
	}

	// S(u, v) = the profile point at u, moved v along the direction
	inline bool TriangulateTrimmedExtrusion(IfcGeometry &geometry, std::vector<IfcBound3D> const& bounds, IfcSurface const& surface, double scaling, int outerTurn = 0)
	{
		auto const& extrusion = surface.ExtrusionSurface;
		if (extrusion.Profile.isComposite || glm::length(extrusion.Direction) == 0) return false;
		auto const& given = extrusion.SweptCurve.points.size() >= 2 ? extrusion.SweptCurve.points : extrusion.Profile.curve.points;
		if (given.size() < 2) return false;
		const glm::dmat4 placement = surface.transformation, back = glm::inverse(placement);
		const glm::dvec3 d = glm::normalize(extrusion.Direction);
		// across the direction the surface is its profile: two axes of that plane
		glm::dvec3 X = glm::cross(d, std::abs(d.x) < 0.9 ? glm::dvec3(1, 0, 0) : glm::dvec3(0, 1, 0));
		X = glm::normalize(X);
		const glm::dvec3 Y = glm::cross(d, X);
		auto profile = std::make_shared<std::vector<glm::dvec3>>(given.begin(), given.end());
		auto line = std::make_shared<TrimPolyline>();
		for (auto const& pt : *profile) line->flat.emplace_back(glm::dot(pt, X), glm::dot(pt, Y));
		double reach = 0;
		for (auto const& pt : *profile) reach = std::max(reach, glm::distance(pt, profile->front()));
		auto curve = std::make_shared<TrimCurve<glm::dvec3>>();
		curve->build(*profile, reach > 0 && glm::distance(profile->front(), profile->back()) < reach * 1e-9);
		auto across = std::make_shared<TrimCurve<glm::dvec2>>();
		across->build(line->flat, reach > 0 && glm::distance(profile->front(), profile->back()) < reach * 1e-9);
		const double last = static_cast<double>(profile->size() - 1);
		const auto local = [=](glm::dvec3 const& pt, double& u, double& v)
		{
			const glm::dvec3 q = glm::dvec3(back * glm::dvec4(pt, 1.0));
			double left = 0;
			u = line->nearest({ glm::dot(q, X), glm::dot(q, Y) }, left);
			u = TrimNearestOnCurve(*across, { glm::dot(q, X), glm::dot(q, Y) }, u, left);
			v = glm::dot(q - curve->at(u), d);
			return left;
		};
		double low = std::numeric_limits<double>::max(), high = std::numeric_limits<double>::lowest();
		for (auto const& bound : bounds)
			for (auto const& pt : bound.curve.points) { double u = 0, v = 0; local(pt, u, v); low = std::min(low, v); high = std::max(high, v); }
		if (!(high > low)) return false;
		TrimSurface s;
		s.u0 = 0; s.u1 = last; s.v0 = low; s.v1 = high;
		double extent = 0;
		for (auto const& pt : line->flat) extent = std::max(extent, glm::distance(pt, line->flat.front()));
		const bool closed = extent > 0 && glm::distance(line->flat.front(), line->flat.back()) < extent * 1e-9;
		s.wrapU = closed;
		s.point = [=](double u, double v)
		{
			double in = u;
			return glm::dvec3(placement * glm::dvec4(curve->at(in) + d * v, 1.0));
		};
		s.project = [=](glm::dvec3 const& pt, glm::dvec2& uv)
		{
			double u = 0, v = 0;
			const double left = local(pt, u, v);
			if (closed) u += std::round((uv.x - u) / last) * last;
			uv = { u, glm::clamp(v, low, high) };
			return left;
		};
		s.direct = true;
		s.outerTurn = outerTurn;
		return TriangulateTrimmedSurface(geometry, bounds, s, scaling);
	}

	// TODO: sent to bimGeometry
    inline void TriangulateBspline(IfcGeometry &geometry, std::vector<IfcBound3D> const& bounds, IfcSurface const& surface, double const scaling, int outerTurn = 0)
	{
		spdlog::debug("[TriangulateBspline({})]");
		Nurbs nurbs{geometry, bounds, surface, scaling};
		nurbs.outerTurn = outerTurn;
		nurbs.fill_geometry();
	}
}