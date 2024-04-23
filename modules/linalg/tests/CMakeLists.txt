# Canonical Durbin's algorithm test
# (Golub & Van Loan, Matrix Computations, 3rd ed)
add_executable(test_durbin_golub_van_loan_3rd_ed
    test_durbin_golub_van_loan_3rd_ed.cpp
    )
target_link_libraries(test_durbin_golub_van_loan_3rd_ed
    INTERFACE b118-linalg-dev
    )
add_test(NAME DurbinGolubVanLoan3rdEdTest
    COMMAND test_durbin_golub_van_loan_3rd_ed
    )

# Canonical Lenvinson's algorithm test
# (Golub & Van Loan, Matrix Computations, 3rd ed)
add_executable(test_levinson_golub_van_loan_3rd_ed
    test_levinson_golub_van_loan_3rd_ed.cpp
    )
target_link_libraries(test_levinson_golub_van_loan_3rd_ed
    INTERFACE b118-linalg-dev
    )
add_test(NAME LevinsonGolubVanLoan3rdEdTest
    COMMAND test_levinson_golub_van_loan_3rd_ed
    )

# Canonical Trench's algorithm test
# (Golub & Van Loan, Matrix Computations, 3rd ed)
add_executable(test_trench_golub_van_loan_3rd_ed
    test_trench_golub_van_loan_3rd_ed.cpp
    )
target_link_libraries(test_levinson_golub_van_loan_3rd_ed
    INTERFACE b118-linalg-dev
    )
add_test(NAME TrenchGolubVanLoan3rdEdTest
    COMMAND test_trench_golub_van_loan_3rd_ed
    )
