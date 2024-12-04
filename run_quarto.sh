#!/bin/bash

# navigate to output directory


quarto render delete_candidates/minimal_quarto_test.qmd 
mv delete_candidates/minimal_quarto_test.html docs/test_test.html