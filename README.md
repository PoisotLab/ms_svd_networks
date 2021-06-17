
## Manuscript versions

[master_draft]: https://poisotlab.github.io/ms_svd_networks/ms_svd_networks-copyedit.pdf
[master_tex]: https://poisotlab.github.io/ms_svd_networks/ms_svd_networks.tex
[master_html]: https://poisotlab.github.io/ms_svd_networks/
[master_preprint]: https://poisotlab.github.io/ms_svd_networks/

<!--- - [:newspaper: preprint][master_preprint] --->
- [:blue_book: website][master_html]
- [:tada: manuscript](https://doi.org/10.3389/fevo.2021.623141)

## A quick primer: 

Generally in ecology we have struggled to define the complexity of ecological networks. Many measures have been proposed, with these usually focusing on the *structural* complexity of the network as opposed to using the ‘*physical*’ complexity of the network. Here we used SVD entropy (based off of information theory) as well as the relative rank deficiency to quantify network complexity. Turns out ecological networks are extremely complex!

Although SVD entropy 'plays well' with the (traditional) structural measures of complexity 
note that entropy values are all above 0.8 where as for the other measures 'complexity 
scores' occupy a larger range between 0 and 1.

![The relationship between SVD entropy and the nestedness (left panel), spectral
radius (central panel) and connectance (right panel) of ecological networks.
Colours indicate the different interaction types of the networks.](figures/others_v_entropy.png){#fig:other}

Pollinator networks are more complex than other types of bipartite networks. Which 
fits in with what what research has been showing us - pollination networks try to 
minimise competition *i.e.* more unique strategies.

![The calculated SVD entropy of different interaction networks of different
interaction types](figures/interactiontype_v_entropy.png){#fig:type}

## About the manuscript building

- Edit `manuscript.md` to make changes to the text.
- The website, draft, and preprint versions are updated after each commit on
  `main`. 
- Previews are also generated for commits on pull-requests, but they need to be
  downloaded as artifacts. The two options to access them are:
  1. Go the `Actions` tab, select the action related to your recent commit, and
     download the `manuscript` zip file in the `Artifacts` section at the bottom
     of the page.
  1. Go the `Checks` tab of your PR (or of one of its commits), select one of
     the jobs on the left, click on `Artifacts` on the right of the page,
     and download the `manuscript` file.

## About the references

- We use [Zotero](https://www.zotero.org/) for references management
- We use the [Better BibTeX](https://retorque.re/zotero-better-bibtex/) plugin
  for citation key generations
- The citation key format we use is meant to convey information on the author,
  date, year, and title. It must be set in the Better BibTeX preferences 
  (`Preferences > Better BibTeX > Citation keys`) as:
    ```
    [auth:fold][year][title:fold:nopunctordash:skipwords:lower:select=1,    1:substring=1,3:capitalize][title:fold:nopunctordash:skipwords:lower:select=2,  2:substring=1,3:capitalize]
    ```
- We also use the same package to automatically export the manuscript references
  to `references.bib`. To do so, right-click on a collection, select Export
  Collection, select the `Better BibTeX` format, and check `Keep updated`.
- It is a good idea to remove a lot of fields that are not strictly speaking
  required for references (`Preferences > Better BibTeX > Export > Fields`). The list of fields we usually ignore is:
    ```
    abstract, copyright, annotation, file, pmid, month, shorttitle, keywords
    ```

## About Markdown formatting

- This is a citation: `@Eckart1936AppOne`. We can also have citations in brackets:
  `[@Eckart1936AppOne]`.
- This is an equation, which we can cite with `@eq:eq1`:
  `$$J = -\frac{1}{\ln(k)}\sum_{i=1}^k s_i\cdot\ln(s_i)$$
  {#eq:eq1}`
- Inline eq. look like this `$J = -\frac{1}{\ln(k)}\sum_{i=1}^k s_i\cdot\ln(s_i)$`
- This is a figure, which we can cite a figure with `@fig:conceptual`:
  `![This is the legend of the figure](figures/conceptual.png){#fig:conceptual}`

More details in the [template repo's `manuscript.md`](https://github.com/PoisotLab/manuscript-template/blob/master/manuscript.md).