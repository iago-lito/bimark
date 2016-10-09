## Current state and RoadMap

This is a first release of the `bimark` package: which will be mostly made of
two parts:

- Superficial "shallow" analysis interface, allowing user to generate,
  manipulate and visualize bilateral mark-recapture data. Written in R.
- Underlying "deep" analysis interface, allowing user to perform a bayesian
  analyis of these data based on their underlying polytope structure. Written in
  R, JAGS and C++.

Althought everything has already been written and run during my past internship,
this first release only ships the "shallow" part. This way, user may:

- assert that..  yes, I am currently working on the overall packaging on this
  code (mostly at nights and weekends).
- try out the first shallow functionnalities:
    - get into the logic
    - provide feedback
    - provide bugreports
    - provide feature requests

In a nutshell, this is an early stage of the project: visit the repo at
[github.repo](github.repo).

## What is Bimark for?

Bimark is dedicated to analysis and simulation of Capture-Mark-Recapture data
in a particular context where the left and right sides of captured
individuals cannot always be matched with one another. We call this the
"bilaterality" problem.

For instance, dolphins may be photo-identified with their dorsal fin. But
depending on whether they exhibit deep marks on the edges or shallow marks on
the side, one may or may not be able to tell whether a left- and a right-picture
correspond to the same fin.

This situation has been formalized in detail by Link *et al.* 2010[^Link2010].
It has been adapted to the problem of bilaterality both by McClintock *et al.*
2013[^McClintock2013] and Bonner *et al.* 2013[^Bonner2013].

The solution we have found to deal with it is mostly inspired from McClintock
and Bonner. Our improved Bayesian sampling algorithm is described in detail in
our M2 report, 2015[^IagNOlivier2015], along with our notations.


[^Link2010]:
    Link *et al.* 2010:
    [doi:10.1111/j.1541-0420.2009.01244.x](https://www.ncbi.nlm.nih.gov/pubmed/19397581)  
[^McClintock2013]:
    McClintock *et al.* 2013:
    [doi:10.1890/12-1613.1](http://onlinelibrary.wiley.com/doi/10.1890/12-1613.1/full)  
[^Bonner2013]:
    Bonner *et al.* 2013:
    [doi:10.1111/biom.12045](http://onlinelibrary.wiley.com/doi/10.1111/biom.12045/abstract)  
[^IagNOlivier2015]:
    soon on [HAL](https://hal.archives-ouvertes.fr/)?

