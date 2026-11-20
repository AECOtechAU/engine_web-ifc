#pragma once

#include <vector>
#include <memory>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <functional>

namespace tinynurbs{
	template<typename T>
	struct RationalSurface;
	typedef RationalSurface<double> RationalSurface3d;
}

namespace webifc::geometry{

	struct IfcGeometry;
	struct IfcBound3D;
	struct BSpline;
	struct IfcSurface;

	constexpr double rotations 		{ 6.0 };
  constexpr auto pi {glm::pi<double>()};
  constexpr auto pi2 {glm::two_pi<double>()};

	// A surface as the trimmed-face mesher needs it
	struct TrimSurface
	{
		double u0 = 0, u1 = 1, v0 = 0, v1 = 1;
		// The surface meets itself across this direction: parameters run on past the end of the range
		bool wrapU = false;
		bool wrapV = false;
		// Look at the surface to find out whether it meets itself
		bool findWraps = false;
		// The nearest point is worked out directly: no starting guess is needed
		bool direct = false;
		// The surface turns by the same amount everywhere, and the bounds lie on it
		bool even = false;
		// The surface is creased wherever this parameter (0: the first, 1: the second) is a whole number
		int creased = -1;
		// The way the outer bound runs in the parameter plane, where the file says: +1 anticlockwise (the face
		// looks the way the surface does), -1 clockwise, 0 not known
		int outerTurn = 0;
		// How far apart two neighbouring directions of the surface may be before another row of points goes in
		double turnStep = 0.35;
		// An inner edge may stray from the surface by this share of its length before it is split
		double stray = 0.066;
		std::function<glm::dvec3(double, double)> point;
		// The parameters of the surface point nearest a position, searched from the parameters given and on
		// the turn of the surface nearest them; returns the distance left
		std::function<double(glm::dvec3 const&, glm::dvec2&)> project;
		glm::dvec2 wrapped(glm::dvec2 const& uv) const;
	};

	// The face inside the bounds, meshed in the surface's parameter plane; false when it could not be
	bool TriangulateTrimmedSurface(IfcGeometry& geometry, std::vector<IfcBound3D> const& bounds, TrimSurface& surface, double scaling);

	class Nurbs{
	
		using uv_point_t = glm::dvec2;
		using uv_points_t = std::vector<uv_point_t>;

		public:

		void fill_geometry();
		explicit Nurbs(IfcGeometry& geometry, std::vector<IfcBound3D>const & bounds, IfcSurface const& surface, double const scaling);
		int outerTurn = 0;

		private:

		void init();
		void fill_geometry_legacy();
		// The trimmed face meshed in the surface's parameter plane; false when it could not be
		bool fill_geometry_trimmed();
		// The parameters of the surface point nearest pt, searched from uv; returns the distance left
		double project(glm::dvec3 const& pt, uv_point_t& uv) const;
		std::vector<double> get_weights() const;
		std::vector<double> get_knots(std::vector<double>const & bs_knots, std::vector<uint32_t> const & bs_mults) const;
		std::vector<glm::dvec3> get_control_points() const;
		uv_points_t get_uv_points() const;
		uv_point_t inverse_evaluation(glm::dvec3 const& pt) const;
		uv_point_t inverse_method(glm::dvec3 const& pt) const;
		std::vector<uint32_t> get_triangulation_uv_points(uv_points_t const& uv_points) const;
		std::vector<double> get_zscores(std::vector<double> const& knots) const;
		std::vector<double> check_knots(std::vector<double> const& knots) const;
		auto get_approximation(glm::dvec3 const& pt, uv_point_t const& range_u, uv_point_t const& range_v) const;

		IfcGeometry& geometry;
		std::vector<IfcBound3D> const& bounds;
		BSpline const& bspline_surface;
		size_t const num_u;
		size_t const num_v;
		double const scaling;
		double minError	{ 0.0001 };
		double maxError	{ 0.01 };
		std::shared_ptr<tinynurbs::RationalSurface3d> nurbs;
		uv_point_t range_knots_u{0.0};
		uv_point_t range_knots_v{0.0};
		glm::dvec3 ptc{0.0};
		glm::dvec3 pth{0.0};
		glm::dvec3 ptv{0.0};
		double dh{0.0};
		double dv{0.0};
		double pr{0.0};
		bool _initialized = false;
		// The surface meets itself across this direction: parameters run on past the end of the range
		bool wrapU = false;
		bool wrapV = false;
		uv_point_t wrapped(uv_point_t const& uv) const;
	};
}