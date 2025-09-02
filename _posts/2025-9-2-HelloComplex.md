---
layout: post
title: Building a Machine Learning System for Protein Complex Prediction
---

Proteins rarely work alone. Life happens when they come together — little machines that fit and unfit, join and dissolve, like a village of parts in constant conversation. Predicting the shape of one protein has been a milestone; predicting how they assemble into complexes feels like the next step. Most fold into shapes that only make sense once you see how they fit together with others.

That’s the project I’ve been working on: building a machine learning system for quaternary structure prediction. Instead of looking at one chain in isolation, the idea is to model how multiple chains interact, forming assemblies that actually do the work of life.

AlphaFold gave us a breathtaking atlas of monomer structures.

There’s no single data stream that tells the whole story. Structures from the PDB, predictions from AlphaFold, interaction maps from cryo-EM, evolutionary signals from sequences — each one is like a partial view through frosted glass. Machine learning, especially graph neural networks, gives us a way to blend these perspectives, to reason across modalities.

The hard part is the scale and fluidity. Complexes can be huge, dynamic, and conditional. You’re always walking a line between overfitting on the static snapshots we have, and missing the living reality where things are constantly in flux.

Still, the effort feels worthwhile. These assemblies are the choreography of biology, and computational models can help us glimpse their dance. It’s equal parts engineering and humility — trying to learn patterns without pretending we can pin them down once and for all.
