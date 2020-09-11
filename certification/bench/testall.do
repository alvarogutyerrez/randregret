clear 
discard 
set more off


cd first 
quietly log using first_test_log, replace
do test
quietly log close
cd ..


cd ..
cd second 
quietly log using second_test_log, replace
do test
quietly log close


cd ..
cd third 
quietly log using third_test_log, replace
do test
quietly log close