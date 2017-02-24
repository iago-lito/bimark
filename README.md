## Current state and RoadMap

This is a first release of the `bimark` package: which will be mostly made of
two parts:

- Superficial "shallow" analysis interface, allowing user to generate,
  manipulate and visualize bilateral mark-recapture data. Written in R.
- Underlying "deep" analysis interface, allowing user to perform a bayesian
  analyis of these data based on their underlying polytope structure. Written in
  R, JAGS and C++.

Althought everything has already been written and run during my past internship
in 2015, this first release only ships the "shallow" part. This way, user may:

- assert that..  yes, I am currently working on the overall packaging on this
  code (mostly at nights and weekends).
- try out the first shallow functionnalities:
    - get into the logic
    - provide feedback
    - provide bugreports
    - provide feature requests

In a nutshell, this is an early stage of the project: visit the repo at
<http://github.com/iago-lito/bimark/>.

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
    M2 [report](http://www.eleves.ens.fr/home/bonnici/Bonnici-Gimenez_2015.pdf)

## Installation

Installing `bimark`, running the tests and building the documentation should be
as easy as:

    > library(devtools)
    > install_github("iago-lito/bimark", build_vignettes=TRUE)

.. but it is [NOT](http://stackoverflow.com/questions/40030414/) XD  
At least not yet. Sorry.

As a temporary solution, please perform this whole development procedure:

    $ git clone https://github.com/iago-lito/bimark
    
Then in R:
    
    > library(devtools)
    > library(testthat)
    > setwd("bimark")
    > document()  
    > use_testthat()
    > setwd("..")
    > install("bimark", dependencies=TRUE, build_vignettes=TRUE)
    > library(bimark)
    > test_package("bimark")
    > run_examples("bimark")

If everything runs fine, then you're done!


## Getting started
    
To generate a model based on simulated data, try:

    > library(bimark)
    > m <- BimarkSimulationModel(N=20, T=5)

The bimark model object is just a list. Access data with `$`:

    > m$n               # number of capture histories actually observed
    > m$LR              # number of observed right-histories
    > m$iOmega          # ids of all histories relevant to these data
    > SeeHist(m$iOmega) # visualize histories and their ids

To generate a model based on actual observation data, use:

    > myData <- example.M
    > m <- BimarkObservationModel(myData)
    > print(m) # unknown number of individuals, since sides haven't been matched
    
To retrieve matrices information from the model, feed dedicated methods with it:

    > GetOmega.B(m) # all unobservable histories that may underlie these data
    > get.A(m)       # observation matrix sorted in polytope order
    > get.B(m)       # kernel of A matrix generating the polytope

Get information with:

    > ?BimarkSimulationModel
    > ?GetOmega.B
    > ?get.A
    
## Troubleshooting

If anything goes wrong (and things *will* go wrong), please file an
issue report on [the repo](https://github.com/iago-lito/bimark/issues).  
We're also pleased to read your feature requests.

## Contributors

Iago-lito (<https://github.com/iago-lito>)  
Dr. Olivier Gimenez, (<olivier.gimenez@cefe.cnrs.fr>)

## License

This package is licensed under the [GPL v3
license](http://www.gnu.org/copyleft/gpl.html). &copy; 2016 Bimark contributors

