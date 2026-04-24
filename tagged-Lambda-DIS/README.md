# tagged-Lambda-DIS
Event generator for leading Lambda baryon tagged DIS process, which is useful for constraining kaon structure function

# Event generation of e p --> e' Lambda X data

**Using ROOT to execute the codes:**
>root -l -q test.cpp

**Using ACLiC to compile and execute:**
>root -l

[0] .x test.cpp+

**Using g++ to compile and execute:**
>sed 's/test()/main()/' test.cpp > test2.cpp

>g++ test2.cpp -o run_sim \`root-config --cflags\` \`root-config --libs\`

>./run_sim


# Reference

Please check G. XIE et al., "Tackling the kaon structure function at EicC",
Chinese Physics C 46 (2022) 6, 064107 [arXiv:2109.08483],
https://inspirehep.net/literature/1923659
for more information about the event generator.


Download ROOT and get started: https://root.cern.ch/

coding issue contact: rwang@impcas.ac.cn

