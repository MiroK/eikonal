#!/bin/bash

ln -T eikonal.x link_to_eikonal.x;
cp link_to_eikonal.x test_results/lin_brent/digits_double;
cp link_to_eikonal.x test_results/lin_brent/digits_double_2;
cp link_to_eikonal.x test_results/lin_brent/digits_double_3;
cp link_to_eikonal.x test_results/lin_geometric;
cp link_to_eikonal.x test_results/lin_newton/digits_double;
cp link_to_eikonal.x test_results/lin_newton/digits_double_2;
cp link_to_eikonal.x test_results/lin_newton/digits_double_3;
rm link_to_eikonal.x


cp compute_rates.py test_results/lin_brent/digits_double;
cp compute_rates.py test_results/lin_brent/digits_double_2;
cp compute_rates.py test_results/lin_brent/digits_double_3;
cp compute_rates.py test_results/lin_geometric;
cp compute_rates.py test_results/lin_newton/digits_double;
cp compute_rates.py test_results/lin_newton/digits_double_2;
cp compute_rates.py test_results/lin_newton/digits_double_3;
