#ifndef INCLUDE_PPJSDM_SATURATED_MODEL
#define INCLUDE_PPJSDM_SATURATED_MODEL

#include <Rcpp.h>

#include "potentials/medium_range_potentials.hpp"
#include "potentials/short_range_potentials.hpp"
#include "../configuration/configuration_manipulation.hpp"
#include "../point/point_manipulation.hpp"
#include "../point/square_distance.hpp"
#include "../utility/for_each_container.hpp"

#include <algorithm> // std::accumulate, std::pop_heap, std::push_heap
#include <functional> // std::greater
#include <memory> // std::shared_ptr
#include <utility> // std::forward, std::move
#include <vector> // std::vector

namespace ppjsdm {

const constexpr char* const short_range_models[] = {
  "exponential",
  "square_exponential",
  "bump",
  "square_bump",
  "Geyer",
  "linear"
};

const constexpr char* const medium_range_models[] = {
  "square_exponential",
  "half_square_exponential",
  "Geyer",
  "linear",
  "half_exponential",
  "exponential",
  "bump",
  "square_bump",
  "tanh"
};

class Saturated_model {
public:
  Saturated_model(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius, unsigned long long int saturation):
  object_(make_short_range_object(model, radius)),
  saturation_(saturation) {}

  Saturated_model(Rcpp::CharacterVector model, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range, unsigned long long int saturation):
  object_(make_medium_range_object(model, medium_range, long_range)),
  saturation_(saturation) {}

  bool is_nonincreasing() const {
    return object_->is_nonincreasing();
  }

  bool is_nonincreasing_after_lower_endpoint() const {
    return object_->is_nonincreasing_after_lower_endpoint();
  }

  bool is_two_valued() const {
    return object_->is_two_valued();
  }

  double apply(double normalized_square_distance, int i, int j) const {
    return object_->apply(normalized_square_distance, i, j);
  }

  double get_square_lower_endpoint(int i, int j) const {
    return object_->get_square_lower_endpoint(i, j);
  }

  unsigned long long int get_saturation() const {
    return saturation_;
  }

  double get_maximum() const {
    return static_cast<double>(saturation_);
  }

private:
  struct Concept {
    virtual ~Concept() {}
    virtual bool is_nonincreasing() const = 0;
    virtual bool is_nonincreasing_after_lower_endpoint() const = 0;
    virtual bool is_two_valued() const = 0;
    virtual double apply(double normalized_square_distance, int i, int j) const = 0;
    virtual double get_square_lower_endpoint(int i, int j) const = 0;
  };

  template<typename Varphi>
  class Concrete_model: public Concept, private Varphi {
  public:
    bool is_nonincreasing() const override {
      return Varphi::is_nonincreasing;
    }

    bool is_nonincreasing_after_lower_endpoint() const override {
      return Varphi::is_nonincreasing_after_lower_endpoint;
    }

    bool is_two_valued() const override {
      return Varphi::is_two_valued;
    }

    template<typename... Args>
    explicit Concrete_model(Args&&... args):
      Varphi(std::forward<Args>(args)...) {}

    double apply(double normalized_square_distance, int i, int j) const override {
      return Varphi::apply(normalized_square_distance, i, j);
    }

    double get_square_lower_endpoint(int i, int j) const override {
      return Varphi::get_square_lower_endpoint(i, j);
    }
  };

  static std::shared_ptr<const Concept> make_short_range_object(Rcpp::CharacterVector model, Rcpp::NumericMatrix radius) {
    const auto model_string(model[0]);
    if(model_string == short_range_models[0]) {
      return std::make_shared<Concrete_model<potentials::Exponential>>(radius);
    } else if(model_string == short_range_models[1]) {
      return std::make_shared<Concrete_model<potentials::Square_exponential>>(radius);
    } else if(model_string == short_range_models[2]) {
      return std::make_shared<Concrete_model<potentials::Bump>>(radius);
    } else if(model_string == short_range_models[3]) {
      return std::make_shared<Concrete_model<potentials::Square_bump>>(radius);
    } else if(model_string == short_range_models[4]) {
      return std::make_shared<Concrete_model<potentials::Strauss>>(radius);
    } else if(model_string == short_range_models[5]) {
      return std::make_shared<Concrete_model<potentials::Linear>>(radius);
    } else {
      Rcpp::stop("Incorrect model entered. A call to show_short_range_models() will show you the available choices.\n");
    }
  }

  static std::shared_ptr<const Concept> make_medium_range_object(Rcpp::CharacterVector model, Rcpp::NumericMatrix medium_range, Rcpp::NumericMatrix long_range) {
    const auto model_string(model[0]);
    if(model_string == medium_range_models[0]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_square_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[1]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_half_square_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[2]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_Geyer>>(medium_range, long_range);
    } else if(model_string == medium_range_models[3]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_linear>>(medium_range, long_range);
    } else if(model_string == medium_range_models[4]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_half_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[5]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_exponential>>(medium_range, long_range);
    } else if(model_string == medium_range_models[6]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_bump>>(medium_range, long_range);
    } else if(model_string == medium_range_models[7]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_square_bump>>(medium_range, long_range);
    } else if(model_string == medium_range_models[8]) {
      return std::make_shared<Concrete_model<potentials::Medium_range_tanh>>(medium_range, long_range);
    } else {
      Rcpp::stop("Incorrect model entered. A call to show_medium_range_models() will show you the available choices.\n");
    }
  }

  std::shared_ptr<const Concept> object_;
  unsigned long long int saturation_;
};

template<typename Point, typename Other>
inline auto apply_potential(const Saturated_model& potential, const Point& point, const Other& other) {
  return potential.apply(normalized_square_distance(point, other), get_type(point), get_type(other));
}

namespace detail {

enum class dispersionMethod {two_values, nonincreasing, nonincreasing_after_lower_endpoint, generic};

template<bool Approximate, dispersionMethod Method>
class compute_dispersion_implementation;

template<bool Approximate>
class compute_dispersion_implementation<Approximate, dispersionMethod::two_values> {
public:
  compute_dispersion_implementation() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&... configurations) const {
    using IntegerType = unsigned long long int;
    using CountType = std::vector<IntegerType>;
    CountType count_vector(number_types);
    CountType deltas(number_types);
    if(static_cast<IntegerType>(varphi.get_saturation()) >= static_cast<IntegerType>(size(configurations...))) {
      for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
        if(!is_equal(current_point, point) && apply_potential(varphi, current_point, point) > 0.) {
          count_vector[get_type(current_point)] += 2;
        }
      }, configurations...);
    } else {
      for_each_container([&deltas, &count_vector, &point, &varphi,
                         saturation = varphi.get_saturation(), &configurations...](const auto& current_point) {
        if(!is_equal(current_point, point)) {
          IntegerType count(0);
          for_each_container([&point, &count, &current_point, &varphi, saturation](const auto& other_point) {
            if(!is_equal(other_point, point) &&
               !is_equal(other_point, current_point) &&
               get_type(other_point) == get_type(point) &&
               count < saturation &&
               apply_potential(varphi, other_point, current_point) > 0.) {
              ++count;
            }
          }, configurations...);
          if(apply_potential(varphi, current_point, point) > 0.) {
            if(count_vector[get_type(current_point)] < saturation) {
              ++count_vector[get_type(current_point)];
            }
            if(count < saturation) {
              ++deltas[get_type(current_point)];
            }
          }
        }
      }, configurations...);
    }
    std::vector<double> dispersion(number_types);
    using size_t = typename decltype(dispersion)::size_type;
    for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
      dispersion[i] = static_cast<double>(count_vector[i]) + static_cast<double>(deltas[i]);
    }
    return dispersion;
  }
};

template<bool Approximate>
class compute_dispersion_implementation<Approximate, dispersionMethod::nonincreasing> {
private:
  using CountType = std::vector<std::vector<double>>;

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(count.size() < varphi.get_saturation()) {
      count.emplace_back(sq);
      std::push_heap(count.begin(), count.end());
    } else if(sq < count[0]) {
      count.emplace_back(sq);
      std::pop_heap(count.begin(), count.end());
      count.pop_back();
    }
  }

  template<typename Vector, typename Point, typename Other>
  static void add_delta(Vector& dispersion, const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(count.size() < varphi.get_saturation()) {
      dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other));
    } else if(sq < count[0]) {
      dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
    }
  }
public:
  compute_dispersion_implementation() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&... configurations) const {
    if(static_cast<unsigned long long int>(varphi.get_saturation())
         >= static_cast<unsigned long long int>(size(configurations...))) {
      CountType count_vector(number_types);
      for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
        if(!is_equal(current_point, point)) {
          compute_dispersion_implementation::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
        }
      }, configurations...);

      std::vector<double> dispersion(number_types);
      using size_t = typename decltype(dispersion)::size_type;
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        dispersion[i] = 2. * std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
          return count + varphi.apply(val, i, get_type(point));
        });
      }
      return dispersion;
    } else {
      CountType count_vector(number_types);
      std::vector<double> dispersion(number_types);

      for_each_container([&dispersion, &count_vector, &point, &varphi, &configurations...](const auto& current_point) {
        if(!is_equal(current_point, point)) {
          std::vector<double> count{};
          for_each_container([&point, &count, &current_point, &varphi](const auto& other_point) {
            if(!is_equal(other_point, point) && !is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
              compute_dispersion_implementation::update_count(varphi, count, current_point, other_point);
            }
          }, configurations...);
          compute_dispersion_implementation::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
          compute_dispersion_implementation::add_delta(dispersion, varphi, count, current_point, point);
        }
      }, configurations...);

      using size_t = typename decltype(dispersion)::size_type;
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        dispersion[i] += std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
          return count + varphi.apply(val, i, get_type(point));
        });
      }
      return dispersion;
    }
  }
};

template<bool Approximate>
class compute_dispersion_implementation<Approximate, dispersionMethod::nonincreasing_after_lower_endpoint> {
private:
  using CountType = std::vector<std::vector<double>>;

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation()) {
        count.emplace_back(sq);
        std::push_heap(count.begin(), count.end());
      } else if(sq < count[0]) {
        count.emplace_back(sq);
        std::pop_heap(count.begin(), count.end());
        count.pop_back();
      }
    }
  }

  template<typename Vector, typename Point, typename Other>
  static void add_delta(Vector& dispersion, const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto sq(normalized_square_distance(point, other));
    if(sq >= varphi.get_square_lower_endpoint(get_type(point), get_type(other))) {
      if(count.size() < varphi.get_saturation()) {
        dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other));
      } else if(sq < count[0]) {
        dispersion[get_type(point)] += varphi.apply(sq, get_type(point), get_type(other)) - varphi.apply(count[0], get_type(point), get_type(other));
      }
    }
  }
public:
  compute_dispersion_implementation() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&... configurations) const {
    if(static_cast<unsigned long long int>(varphi.get_saturation()) >= static_cast<unsigned long long int>(size(configurations...))) {
      CountType count_vector(number_types);
      for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
        if(!is_equal(current_point, point)) {
          compute_dispersion_implementation::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
        }
      }, configurations...);

      std::vector<double> dispersion(number_types);
      using size_t = typename decltype(dispersion)::size_type;
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        dispersion[i] = 2. * std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
          return count + varphi.apply(val, i, get_type(point));
        });
      }
      return dispersion;
    } else {
      CountType count_vector(number_types);
      std::vector<double> dispersion(number_types);

      for_each_container([&count_vector, &dispersion, &point, &varphi, configurations...](const auto& current_point) {
        if(!is_equal(current_point, point)) {
          std::vector<double> count{};
          for_each_container([&point, &count, &current_point, &varphi](const auto& other_point) {
            if(!is_equal(other_point, point) && !is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
              compute_dispersion_implementation::update_count(varphi, count, current_point, other_point);
            }
          }, configurations...);
          compute_dispersion_implementation::update_count(varphi, count_vector[get_type(current_point)], current_point, point);
          compute_dispersion_implementation::add_delta(dispersion, varphi, count, current_point, point);
        }
      }, configurations...);

      using size_t = typename decltype(dispersion)::size_type;
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        dispersion[i] += std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0., [&varphi, i, &point](double count, const auto& val) {
          return count + varphi.apply(val, i, get_type(point));
        });
      }
      return dispersion;
    }
  }
};

template<bool Approximate>
class compute_dispersion_implementation<Approximate, dispersionMethod::generic> {
private:
  using CountType = std::vector<std::vector<double>>;

  template<typename Point, typename Other>
  static void update_count(const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    if(count.size() < varphi.get_saturation()) {
      count.emplace_back(disp);
      std::push_heap(count.begin(), count.end(), std::greater<double>{});
    } else if(disp > count[0]) {
      count.emplace_back(disp);
      std::pop_heap(count.begin(), count.end(), std::greater<double>{});
      count.pop_back();
    }
  }

  template<typename Vector, typename Point, typename Other>
  static void add_delta(Vector& dispersion, const Saturated_model& varphi, typename CountType::value_type& count, const Point& point, const Other& other) {
    const auto disp(apply_potential(varphi, other, point));
    if(count.size() < varphi.get_saturation()) {
      dispersion[get_type(point)] += disp;
    } else if(disp > count[0]) {
      dispersion[get_type(point)] += disp - count[0];
    }
  }
public:
  compute_dispersion_implementation() {}

  template<typename Point, typename... Configurations>
  auto operator()(const Saturated_model& varphi,
                const Point& point,
                R_xlen_t number_types,
                Configurations&... configurations) const {
    if(static_cast<unsigned long long int>(varphi.get_saturation()) >= static_cast<unsigned long long int>(size(configurations...))) {
      CountType count_vector(number_types);
      for_each_container([&count_vector, &point, &varphi](const auto& current_point) {
        if(!is_equal(current_point, point)) {
          compute_dispersion_implementation::update_count(varphi, count_vector[get_type(current_point)], point, current_point);
        }
      }, configurations...);

      std::vector<double> dispersion(number_types);
      using size_t = typename decltype(dispersion)::size_type;
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        dispersion[i] = 2. * std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0.);
      }
      return dispersion;
    } else {
      CountType count_vector(number_types);
      std::vector<double> dispersion(number_types);

      for_each_container([&count_vector, &dispersion, &point, &varphi, configurations...](const auto& current_point) {
        if(!is_equal(current_point, point)) {
          std::vector<double> count{};
          for_each_container([&point, &count, &current_point, &varphi](const auto& other_point) {
            if(!is_equal(other_point, point) && !is_equal(other_point, current_point) && get_type(other_point) == get_type(point)) {
              compute_dispersion_implementation::update_count(varphi, count, other_point, current_point);
            }
          }, configurations...);
          compute_dispersion_implementation::update_count(varphi, count_vector[get_type(current_point)], point, current_point);
          compute_dispersion_implementation::add_delta(dispersion, varphi, count, current_point, point);
        }
      }, configurations...);

      using size_t = typename decltype(dispersion)::size_type;
      for(size_t i(0); i < static_cast<size_t>(number_types); ++i) {
        dispersion[i] += std::accumulate(count_vector[i].begin(), count_vector[i].end(), 0.);
      }
      return dispersion;
    }
  }
};

} // namespace detail

template<bool Approximate = false, typename Point, typename... Configurations>
inline auto compute_dispersion(const Saturated_model& model,
                               const Point& point,
                               R_xlen_t number_types,
                               Configurations&&... configurations) {
  if(model.is_two_valued()) {
    return detail::compute_dispersion_implementation<Approximate, detail::dispersionMethod::two_values>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing()) {
    return detail::compute_dispersion_implementation<Approximate, detail::dispersionMethod::nonincreasing>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else if(model.is_nonincreasing_after_lower_endpoint()) {
    return detail::compute_dispersion_implementation<Approximate, detail::dispersionMethod::nonincreasing_after_lower_endpoint>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  } else {
    return detail::compute_dispersion_implementation<Approximate, detail::dispersionMethod::generic>{}(model, point, number_types, std::forward<Configurations>(configurations)...);
  }
}

} // namespace ppjsdm

#endif // INCLUDE_PPJSDM_SATURATED_MODEL
