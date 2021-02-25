# Reviewer 1

## Please summarize the main findings of the study.

Strydom and Poisot make an interesting analysis about the complexity of the
ecological network. In particular, they study the interaction between entropy,
as derived from the Shannon index and other variables, like richness,
interaction type, etc. One of their conclusion is that the singular value
decomposition entropy can be used for unifying the concepts of complexity in the
ecological networks. They find that the complexity of the ecological networks is
extremely large.

## Please highlight the limitations and strengths.

The paper is presented in a clear manner, it has an excellent introduction and
bibliographic literature. Figures are well-drawn and clear.

I see two main limitations. The first is that it is not defined a specific
hypothesis, and all the paper stands in between a theoretical paper and an
exercise. The second is that the source of data comes from a single repository
organized by the same author of the software they use for computation. We,
therefore, must somewhat rely upon these data, without any chance to test the
concepts expressed with independent datasets. As a reader, I am therefore unable
to assess to which extent the conclusions are supported by the data.

> We have expanded on the data and methods section to give a more detailed
> overview of the bipartite networks used, which we present as a summary table.
> This includes the sample size, latitudinal range, and top and bottom richness
> for the different interaction types to highlights the diversity of datasets
> used

## Please comment on the methods, results and data interpretation. If there are any objective errors, or if the conclusions are not supported, you should detail your concerns.

As expressed above, I would recommend more clarity about the used dataset, and
possibly I would add data from other sources. Alternatively, I would recommend
limiting the analysis to theoretical considerations.

It is not clear which of the data present in the cited repository have been
used, and who exactly prepared the database, when, and with which observational
method.

## Other comments

I found a single grammar mistake, repeated twice: 'defficiency' instead of
deficiency.

> This has been corrected

# Reviewer 2

## Please summarize the main findings of the study.

Remarkable is, in my opinion, the result that complexity seems to be negatively
correlated to nestedness, to the spectral radius, and connectance, highlighting
that mutualistic networks tend to be more complex in order to minimize
competition relationships by reducing the overlap of interactions. Therefore,
the lack of complexity of an ecological network can be measured both through
nestedness and connectance.

## Please highlight the limitations and strengths.

To make the text more understandable even for non-specialized readers, I suggest
the authors give some important information in the appendix. Authors should
better formalise and clarify:

a) the concept of simulated annealing which, although widely used in various
fields, may not be known and above all a greater detail of the steps necessary
to generate networks with the highest, or lowest, possible SVD entropy values.

b) some details on the nomenclature of the matrices used (for example singular
matrix, diagonal matrix, transposed matrix, etc.) so that the reader has a
clearer basis for calculating the entropy VSD.

> This has been addressed in-text as opposed to an appendix. Although Singular
> Value Decomposition forms the initial step in calculating entropy it is the
> Sigma matrix that is of interest for further calculations and it is more
> important that the reader grasp those concepts as opposed to getting bogged
> down in some of the other mathematical details. That being said we have
> briefly included how the U and VT matrices are derived (ca. l. 152 in .md) and
> expanded on the definition of a diagonal matrix (ca. l. 153 in .md)

c) how is nestedness calculated?

> We already point to the paper which describes calculating nestedness in l. 213
> but have also added that “nestedness is calculated based on the number of
> shared links between species pairs...” in the preceding paragraph (ca l. 219)
> to better mirror how we have described connectance and spectral radius.

## Please comment on the methods, results and data interpretation. If there are any objective errors, or if the conclusions are not supported, you should detail your concerns.

The manuscript is well written and clear in the discussion of the results.

## Please provide your detailed review report to the editor and authors (including any comments on the Q4 Check List):

The manuscript is focused on the application of Singular Value Decomposition
(SVD) Entropy to 220 bipartite networks of different types of interactions for
the definition of complexity in relation to the consideration of the "physical"
complexity of the network rather than on the "behavioural" complexity. The
subject is of great interest to the ecological sciences and this manuscript
offers a series of results and considerations that undoubtedly lead to a serious
advancement of knowledge in this area.

However, to make the text more understandable even for non-specialized readers,
I suggest the authors give some important information in the appendix. Authors
should better formalise and clarify:

a) the concept of simulated annealing which, although widely used in various
fields, may not be known and above all a greater detail of the steps necessary
to generate networks with the highest, or lowest, possible SVD entropy values.

b) some details on the nomenclature of the matrices used (for example singular
matrix, diagonal matrix, transposed matrix, etc.) so that the reader has a
clearer basis for calculating the entropy VSD.

c) how is nestedness calculated?

