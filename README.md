# NeuralProgenitors
<p>
  This repository contains the scripts used to generate the data presented in
  <em>“Neural Progenitors as a Novel Pathogenic Mechanism in Microcephaly.”</em>
</p>

<p>The following scripts:</p>
<ul>
  <li><code>cortical_analysis.R</code></li>
  <li><code>rostal_analysis.R</code></li>
</ul>
<p>were used to analyze the scRNA-seq data using Seurat.</p>

<p>
  The script <code>read_human_embryonic_atlas.R</code> was used to process the
  data from Braun et al. (2023), <em>Science</em> (downloaded from
  <a href="https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077">
    https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077
  </a>),
  and to generate the scRNA reference matrix for CIBERSORT bulk deconvolution.
</p>

<p>
  The script <code>label_transfer.R</code> was used to annotate the cortical
  scRNA-seq dataset using the Braun et al. dataset as a reference.
</p>
`
