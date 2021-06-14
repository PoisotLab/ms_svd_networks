
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