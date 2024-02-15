/*   libb118
 *
 *   modules/frlap/src/gdm/calculators/simple.cpp
 *
 *   Simple fractional Laplacian calculator using general differences
 *   approximations
 *
 *   Copyright (C) 2023  Guilherme F. Fornel <gffrnl@gmail.com>
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <b118/frlap/gdm/calculators/simple.hpp>
#include <b118/closest_sorted.hpp>
#include <b118/linspace.hpp>

using b118::frlap::gdm::calculators::Problem;

class Outputter;

///
/// @brief Function template to compute the normalization constant
///        associated to the singular integral representation of
///        the fractional Laplacian.
///
/// @param probptr     Pointer to the problem.
/// @param fnameprefix Problem's filenames prefix
/// @param fpathprefix Problem's filepaths prefix
///
/// @return Maximum absolute error
///
double compute_frlap_for_problem(std::unique_ptr<Problem> & probptr    ,
                                 std::string                fnameprefix,
                                 std::string                fpathprefix) {
    double maxabserr = 0;

    std::ofstream ofs_info;  // fstream to output the problem's information
    std::ofstream ofs_dat;   // fstream to output the problem's results

    std::string const ofname_info = fnameprefix + ".info";
    std::string const ofname_dat  = fnameprefix + ".dat";
    std::string const ofpath_info = fpathprefix + ".info";
    std::string const ofpath_dat  = fpathprefix + ".dat";

    // Problem parameters:
    double      const ealpha = probptr->ealpha;
    double      const a      = probptr->domain.a;
    double      const b      = probptr->domain.b;
    std::size_t const numnod = probptr->numnod;
    double      const xa     = probptr->viswin.a;
    double      const xb     = probptr->viswin.b;

    // Grid parameters:
    double deltax;  // grid spacing
    std::size_t ja, jb;
    std::size_t numnod_xa_xb;

    // Solution Vectors:
    std::vector<double> X;  // grid points
    std::vector<double> Y;  // grid function values
    std::vector<double> FLY0;
    std::vector<double> FF;

    std::cout << "\n\t computing the fractional Laplacian for problem "
              << fnameprefix
              << " ...\n"
              << std::endl;

    // Creation of the grid values and grid function values
    X = b118::linspace<std::vector, double>(a, b, numnod);
    Y.resize(X.size());
    // std::transform(std::execution::par,
    //                X.cbegin(), X.cend(),
    //                Y.begin(),
    //                [&probptr](double x) -> double { return probptr->y(x); });
    std::transform(X.cbegin(), X.cend(),
                   Y.begin(),
                   [&probptr](double x) -> double { return probptr->y(x); });

    ja = b118::closest_sorted(X.cbegin(), X.cend(), xa);  // throwns exception
    jb = b118::closest_sorted(X.cbegin(), X.cend(), xb);  // throwns exception

    // Print domain parameters to the console
    (void) std::puts(" Domain parameters:");
    (void) std::printf("\t[a, b] = [%+.2f, %+.2f]\n", a, b);
    (void) std::printf("\tn = %ld\n", numnod);
    (void) std::printf("\th = %.5f\n", deltax);
    (void) std::puts("\t");

    // Actual computation
    probptr->method_ptr->ealpha = ealpha;
    probptr->method_ptr->deltax = deltax;
    probptr->method_ptr->compute(Y, ja, jb, FLY0);
    numnod_xa_xb = FLY0.size();
    FF.resize(numnod_xa_xb);
    probptr->compute_far_field(X, ja, jb, FF.data());


    {
        // Write the results to a file:
        ofs_dat.open(ofpath_dat);
        if (!ofs_dat.is_open()) {
            std::cerr << "error: could not open file "
                      << ofpath_dat
                      << " for writing."
                      << std::endl;
            throw "file writting error";
        }

        ofs_dat << std::left << std::setw(24) << "X";
        ofs_dat << std::left << std::setw(24) << "Y";
        ofs_dat << std::left << std::setw(24) << "FLY0";
        ofs_dat << std::left << std::setw(24) << "FF";
        ofs_dat << std::left << std::setw(24) << "FLY";
        ofs_dat << std::left << std::setw(24) << "Exact";
        ofs_dat << std::left << std::setw(24) << "AbsError";
        ofs_dat << '\n';

        for (size_t j = 0; j < numnod_xa_xb; ++j) {
            double approx, exact;

            ofs_dat << std::scientific << std::setprecision(14) << std::showpos;
            ofs_dat << std::left << std::setw(24) << X[j+ja];
            ofs_dat << std::left << std::setw(24) << Y[j+ja];
            ofs_dat << std::left << std::setw(24) << FLY0[j];
            ofs_dat << std::left << std::setw(24) << FF[j];

            approx = FLY0[j]+FF[j];
            exact  = probptr->frLap_y(X[j+ja]);

            ofs_dat << std::left << std::setw(24) << approx;
            ofs_dat << std::left << std::setw(24) << exact;

            {
                double const abserr = std::fabs(exact-approx);
                ofs_dat << std::left << std::setw(24) << abserr;
                if (abserr > maxabserr) maxabserr = abserr;
            }

            ofs_dat << std::endl;
        }

        ofs_dat.close();
    }

    {
        // Problem's information file
        ofs_info.open(ofpath_info);
        if (!ofs_info.is_open()) {
            std::cerr << "error: could not open file "
                      << ofpath_info
                      << " for writing."
                      << std::endl;
            throw "file writting error";
        }
        ofs_info << "# problem " << fnameprefix << " information file\n"
                 << std::endl;
        ofs_info << "output filename: " << ofname_dat
                 << std::endl;
        ofs_info << "fractional Laplacian exponent: " << ealpha
                 << std::endl;
        ofs_info << std::scientific << std::setprecision(5);
        ofs_info << "domain: [" << a << ", " << b << "]"
                 << std::endl;
        ofs_info << "vis. window: [" << xa << ", " << xb << "]"
                 << std::endl;
        ofs_info << "number of nodes in domain: " << numnod
                 << std::endl;
        ofs_info << "number of nodes in vis. window: " << numnod_xa_xb
                 << std::endl;
        ofs_info << "grid spacing: " << std::setprecision(14) << deltax
                 << std::endl;
        ofs_info.close();
    }

    return maxabserr;
}


class Outputter {
    class OFileExec {
        char const * const suffix;

     protected:
        std::ofstream fs;
        std::string name;
        std::string path;

     public:
        OFileExec(char const * suffix, std::string const & execname,
                                       std::string const & prefix)
            : suffix(suffix),
              name(execname + suffix),
              path(prefix + suffix)
        {}

        virtual ~OFileExec() {}

        std::string filename() const { return name; }

        void open() {
            fs.open(path);
            if (!fs.is_open()) {
                std::cerr << "error: could not open file "
                          << path
                          << " for writing."
                          << std::endl;
                throw "could not open";
            }
        }

        void close() { fs.close(); }

        bool is_open() const { return fs.is_open(); }


        friend
        OFileExec& operator<<(OFileExec& ofile, char const & c) {
            ofile.fs << c;
            return ofile;
        }

        friend
        OFileExec& operator<<(OFileExec& ofile, char const * str) {
            ofile.fs << str;
            return ofile;
        }

        friend
        OFileExec& operator<<(OFileExec& ofile, std::string const & str) {
            ofile.fs << str;
            return ofile;
        }
    };

    class OFileExecInfo final : public OFileExec {
     public:
        explicit OFileExecInfo(std::string const & execname,
                               std::string const & prefix)
            : OFileExec(".info", execname, prefix)
        {}

        OFileExecInfo() : OFileExecInfo("", "") {}

        OFileExecInfo& operator=(OFileExecInfo&& other) {
            fs   = std::move(other.fs);
            name = std::move(other.name);
            path = std::move(other.path);
            return *this;
        }
    };

    class OFileExecMaxAbsErr final : public OFileExec {
     public:
        explicit OFileExecMaxAbsErr(std::string const & execname,
                                    std::string const & prefix)
            : OFileExec(".info", execname, prefix)
        {}

        OFileExecMaxAbsErr() : OFileExecMaxAbsErr("", "") {}

        OFileExecMaxAbsErr& operator=(OFileExecMaxAbsErr&& other) {
            fs   = std::move(other.fs);
            name = std::move(other.name);
            path = std::move(other.path);
            return *this;
        }
    };

    struct ProbFile {
        std::string nameprefix;
        std::string pathprefix;
    };

 public:
    OFileExecInfo info;
    OFileExecMaxAbsErr maxabserr;
    std::vector<ProbFile> prob;

    explicit Outputter(std::string const & execname,
     std::vector<std::unique_ptr<Problem>> const & prob_ptrs) {
        std::string const prefix  = get_prefix(execname);
        std::size_t const numprob = prob_ptrs.size();

        info = OFileExecInfo(execname, prefix);
        maxabserr = OFileExecMaxAbsErr(execname, prefix);

        prob.resize(numprob);
        for (std::size_t p = 0; p < numprob; ++p)
            prob[p].nameprefix = execname + "_" + (prob_ptrs[p]->label);
        for (std::size_t p = 0; p < numprob; ++p)
            prob[p].pathprefix = prefix   + "_" + (prob_ptrs[p]->label);
    }

    std::size_t num_problems() const { return prob.size(); }

 protected:
    // TODO(gffrnl): Check if it works on Windows
    static std::string get_prefix(std::string const & execname) {
        // Only works on C++17 or newer (__cplusplus >= 201703L)
        namespace fs = std::filesystem;

        // Try to create a folder with the executable name in the current path
        std::string outprefix{""};
#ifdef _WIN32
        char const * c = "\\";
#else
        char const * c = "/";
#endif
        std::string dirpath;
        dirpath = fs::current_path().string()
                                    .append(c)
                                    .append(execname)
                                    .append("_out");
        fs::create_directories(dirpath);

        // If the folder was created, then append its name
        if (fs::exists(dirpath))
            outprefix.append(dirpath).append(c);

        outprefix.append(execname);
        return outprefix;
    }
};

std::string get_execname(char const * argv0);

int main(int argc, char * argv[]) {
    using b118::frlap::gdm::calculators::prob_ptrs;

    std::string const execname = get_execname(argv[0]);  // The executable name.

    Outputter out(execname, prob_ptrs);

    std::cout << "\n  *** frlap_gdm_calculator ***\n\n"
              << "  number of problems: "
              << out.num_problems()
              << std::endl;

    out.info.open();
    out.info << "# " << execname << " information file\n\n";
    out.info << "maximum abs. error file: "
             << out.maxabserr.filename()
             << '\n';
    out.info << "number of problems: " << out.num_problems() << '\n';
    out.info << "problems:\n";

    out.maxabserr.open();
    out.maxabserr << "problem, maxabserr\n";
    for (std::size_t p = 0; p < out.num_problems(); ++p) {
        double maxabserr;
        maxabserr =
            compute_frlap_for_problem(prob_ptrs[p],
                                      out.prob[p].nameprefix,
                                      out.prob[p].pathprefix);

        out.maxabserr << out.prob[p].nameprefix
                      << ", "
//                      << std::scientific
//                      << std::setprecision(14)
                      << maxabserr
                      << 'n';
        out.info << "\t" << out.prob[p].nameprefix << '\n';
    }
    out.maxabserr.close();


    out.info << "\nprogram termination with success\n";
    out.info.close();

    return EXIT_SUCCESS;
}


//
// @brief Gets the current executable name
//
// TODO(gffrnl): Check if it works on Windows
//
// @param argv0 The C string in `argv[0]` from `main(int argc, char * argv)`
//
// @return A std::string with the executable name
//
std::string get_execname(char const * argv0) {
#ifdef _WIN32
    char const * c = "\\";
#else
    char const * c = "/";
#endif
    std::string s(argv0);

    if (s.length() == 0)
        throw std::invalid_argument("get_execname(): "
                                    "argument is an empty C string");

    std::string::size_type const n = std::string(argv0).rfind(c);

    if (n + 1 == s.length())
#ifdef _WIN32
        throw std::invalid_argument("get_execname(): "
                                    "argument `\\` terminated");
#else
        throw std::invalid_argument("get_execname(): "
                                    "argument `/` terminated");
#endif

  s = s.substr(n + 1);

  return s.substr(0, s.rfind('.'));
}
