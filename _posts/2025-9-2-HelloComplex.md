---
layout: post
title: Mapping the Quaternary Landscape
---
### Building a Machine Learning System for Protein Complex Prediction

Proteins rarely work alone. Life happens when they come together — little machines that fit and unfit, join and dissolve, like a village of parts in constant conversation. Predicting the shape of one protein has been a milestone, AlphaFold gave us a breathtaking atlas of monomer structures; predicting how they assemble into complexes feels like the next step. Most fold into shapes that only make sense once you see how they fit together with others.

That’s the project the team at Jäntra has been working on: building a machine learning system for quaternary structure prediction. Instead of looking at one structure in isolation, the idea is to model how multiple structures interact, forming assemblies that actually do the work of life.

![Dalton PDB](/images/dalton_beta.gif "Editing Molecular Structure")

There’s no single data stream that tells the whole story. Structures from the PDB, predictions from AlphaFold, interaction maps from cryo-EM, evolutionary signals from sequences — each one is like a partial view through frosted glass. Machine learning, especially deep neural networks, gives us a way to blend these perspectives, to reason across modalities. The hard part is the scale and fluidity. Complexes can be huge, dynamic, and conditional. 

Our starting approach and can be summerized by the following bottom-up strategy for quaternary structure prediction: 

1. Start with a closed world → pick an organism like E. coli, where we already have tons of experimental + predicted structures in the AlphaFold DB.

2. Bucketize monomers → cluster proteins into functional or structural groups (could be sequence families, Pfam domains, structural folds, GO annotations, etc.). This helps tame the combinatorial explosion — you’re not checking every protein against every other.

3. Test for intra-bucket complexes → proteins in the same bucket are more likely to form homomers or closely related assemblies. You can start small here, validating known complexes (ribosomal proteins, polymerases, etc.).

4. Expand to inter-bucket pairs → once you have confidence, test cross-bucket interactions. This is where you might discover novel or less obvious complexes.

5. Layer in ML → the buckets become priors or constraints for your machine learning model: instead of feeding the model the entire proteome, you give it “likely interaction neighborhoods.”

This approach we believe is similar to reducing a huge search space into manageable islands, then letting the model explore the coastlines where those islands might connect.

It's a very different framing from AlphaFold-Multimer (which is trained end-to-end on complexes), but it’s a neat exploratory approach that could complement it: let the data itself suggest which assemblies to test.

You’re always walking a line between overfitting on the static snapshots we have, and missing the living reality where things are constantly in flux. Still, the effort feels worthwhile. These assemblies are the choreography of biology, and computational models can help us glimpse their dance. It’s equal parts engineering and humility — trying to learn patterns without pretending we can pin them down once and for all.
