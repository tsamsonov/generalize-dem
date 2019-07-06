#include<vector>
#include<cmath>
#include<iostream>
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>

namespace py = pybind11;

std::vector<unsigned int> find_up_cell(double *inr, unsigned int nrow, unsigned int ncol, unsigned int i, unsigned int j) {

    std::vector<std::vector<double>> w = {
        {0.70710678, 1, 0.70710678},
        {1, 1, 1},
        {0.70710678, 1, 0.70710678}
    };

    std::vector<int> shift = {-1, 0, 1};

    double minmax = 4e9;
    auto kmin = 1;
    auto lmin = 1;

    auto k1 = (i > 0) ? 0 : 1,
         k2 = (i < nrow-1) ? 3 : 2,
         l1 = (j > 0) ? 0 : 1,
         l2 = (j < ncol-1) ? 3 : 2;

    for (int k = k1; k < k2; k++) {
        auto ik = i + shift[k];
        for (int l = l1; l < l2; l++) {
            auto jl = j + shift[l];
            auto temp = (inr[i * ncol + j] - inr[ik * ncol + jl]) * w[k][l];

            if ((temp > 0) && (temp < minmax)) {
                minmax = temp;
                kmin = k;
                lmin = l;
            }
        }
    }

    std::vector<unsigned int> upcell = {i + shift[kmin], j + shift[lmin]};

    return upcell;
}

py::array_t<double> extract_streams(py::array_t<double> in_raster, py::array_t<double> out_raster, double min_acc, int min_len) {

    auto buf_in = in_raster.request(), buf_out = out_raster.request();

    auto nrow = buf_in.shape[0];
    auto ncol = buf_in.shape[1];

    auto *inptr = (double *) buf_in.ptr,
         *outptr = (double *) buf_out.ptr;

    std::vector<std::vector<unsigned int>> stream(min_len);

    for (unsigned int i = 0; i < nrow; i++) {
        std::cout << "Processing line " << i << " from " << nrow << std::endl;
        for (unsigned int j = 0; j < ncol; j++) {

            if (inptr[i * ncol + j] < min_acc)
                continue;

            auto ik = i, jk = j;
            auto n = 0;
            auto stream_reached = false;

            auto idx = 0;

            while (n < min_len) {
                stream[n] = std::vector<unsigned int>{ik, jk};

                auto upcell = find_up_cell(inptr, nrow, ncol, ik, jk);
                auto ik_new = upcell[0], jk_new = upcell[1];

                if ((ik_new == ik) && (jk_new == jk))
                    break;

                idx = ik_new * ncol + jk_new;

                if (inptr[idx] < min_acc)
                    break;

                if (outptr[idx] > 0) { // stream reached
                    stream_reached = true;
                    n++;
                    break;
                }

                ik = ik_new;
                jk = jk_new;
                n++;

            }


            if (n == min_len || stream_reached) {
                for (int k = 0; k < n; k++) {
                    outptr[stream[k][0] * ncol + stream[k][1]] = 1;
                }

                if (!stream_reached) {
                    outptr[idx] = 1;
                    while (true) {

                        auto upcell = find_up_cell(inptr, nrow, ncol, ik, jk);
                        auto ik_new = upcell[0], jk_new = upcell[1];

                        if ((ik_new == ik) && (jk_new == jk))
                            break;

                        idx = ik_new * ncol + jk_new;

                        if (inptr[idx] < min_acc)
                            break;

                        if (outptr[idx] > 0) { // stream reached
                            break;
                        }

                        outptr[idx] = 1;

                        ik = ik_new;
                        jk = jk_new;
                    }
                }
            }
        }
    }

    return out_raster;
}

PYBIND11_MODULE(StreamExtractor3, m) {
    m.doc() = R"pbdoc(
        C++ plugin for extracting streams from raster DEM
        -----------------------

        .. currentmodule:: StreamExtractor

        .. autosummary::
           :toctree: _generate

           extract_streams

    )pbdoc";

    m.def("extract_streams", &extract_streams, R"pbdoc(
        Extract streams from flow accumulation raster using minimum length and minimum flow accumulation criteria
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}