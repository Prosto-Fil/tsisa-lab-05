#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using std::cout;
using std::endl;

class neuron {
private:
    double x;
    double y;

public:
    neuron() : x(0), y(0) {};

    void set_x(double buf_x) {
        x = buf_x;
    }

    void set_y(double buf_y) {
        y = buf_y;
    }

    auto get_x() const noexcept -> double {
        return x;
    }

    auto get_y() const noexcept -> double {
        return y;
    }
};

auto random(double min, double max) -> double {
    return (double)(rand()) / RAND_MAX * (max - min) + min;
}

auto linear_function(double c, double d, double x) -> double {
    return  c * x + d;
}

auto random_error(double a, double b, double A) -> double {
    return A * random(a, b);
}

auto square_func(const std::vector<neuron>& neurons, double c, double d, double N) -> double {
    double sum = 0;
    for (auto i = 0; i < N; i++) {
        auto x = neurons[i].get_x();
        auto y = neurons[i].get_y();
        double new_y = linear_function(c, d, x);
        sum += pow(new_y - y, 2);
    }
    return sum;
}

auto fill_network(double c, double d, double A, double a, double b, size_t neurons_number) -> std::vector<neuron> {
    auto neuron_distance = (b - a) / static_cast<double>(neurons_number + 1);
    std::vector<neuron> neurons(neurons_number);
    auto current_x = a;
    for (auto& neuron : neurons) {
        current_x += neuron_distance;
        neuron.set_x(current_x);
        neuron.set_y(linear_function(c, d, neuron.get_x()) + random_error(-0.5, 0.5, A));
    }
    return neurons;
}

auto passive_search(const std::vector<neuron>& neurons, double a, double b, double N) -> double {
    auto min_d = neurons[0].get_y();
    auto max_d = neurons[0].get_y();
    for (const auto& neuron : neurons) {
        if (neuron.get_y() < min_d) {
            min_d = neuron.get_y();
        }
        if (neuron.get_y() > max_d) {
            max_d = neuron.get_y();
        }
    }
    auto d = random(min_d, max_d);

    auto max_c = 1.26;
    auto min_c = -1.37;
    std::vector<double> vec_c;
    const int number_of_iterations = 99;
    std::vector<double> sum;
    for (size_t k = 0; k < number_of_iterations; k++) {
        vec_c.push_back(min_c + ((max_c - min_c) / static_cast<double>(number_of_iterations + 1)) * (k + 1));
        sum.push_back(square_func(neurons, vec_c[k], d, N));
    }

    auto min_sum = std::min_element(sum.begin(), sum.end());
    auto num = std::distance(sum.begin(), min_sum);
    return vec_c[num];
}

auto dihotomia(const std::vector<neuron>& neurons, double a, double b, double c, double N) -> double {
    const double eps = 0.1;
    const double delta = 0.01;
    do {
        auto d_left = 0.5 * (b + a) - delta;
        auto d_right = 0.5 * (b + a) + delta;
        auto f_left = square_func(neurons, c, d_left, N);
        auto f_right = square_func(neurons, c, d_right, N);
        if (f_right > f_left) {
            b = d_right;
        } else {
            a = d_left;
        }
    } while ((b - a) > eps);
    return ((b + a) / 2);
}

void print(const std::vector<neuron>& neurons) {
    cout << "-------------------" << endl;
    cout << "| " << std::setw(5) << std::left << "x" << " | " << std::setw(8) << std::left << "y" << "|" << endl;
    cout << "-------------------" << endl;
    for (const auto& neuron : neurons) {
        cout << "| " << std::setw(5) << std::left << neuron.get_x() << " | " << std::setw(8) << std::left << neuron.get_y() << "|" << endl;
    }
    cout << "-------------------" << endl;
}

int main() {
    const double a = -4;
    const double b = 2;
    const double c = 0;
    const double d = 3;
    const size_t N = 24;
    const double A = 0.1;

    auto network = fill_network(c, d, A, a, b, N);

    auto c_find = passive_search(network, a, b, N);

    cout << "c = " << c_find << endl;
    cout << "d = " << dihotomia(network, a, 4, c_find, N) << endl;

    cout << endl << "Network: " << endl;

    print(network);

    return 0;
}
