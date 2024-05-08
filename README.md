EDITS IN PROGRESS....

**Background** 

Precision functional mapping (PFM) refers to a suite of approaches for delineating functional brain areas and networks at the individual level.

**Example PFM dataset**  

An example PFM dataset (“ME01”, study author CJL) is available online at OpenNeuro.org (https://openneuro.org/datasets/ds005118/) in part for the purpose of providing example data for this tutorial. Please cite Lynch et al., 2024 (Citation TBD) when using this code, and Lynch et al. 2020 Cell Reports (https://pubmed.ncbi.nlm.nih.gov/33357444/) when publishing work that includes the ME01 dataset.

**Usage**

The goal of this tutorial is to describe step-by-step how PFM is performed in our lab. This tutorial focuses on the ME01 dataset, but in principle the code can be adapted to accommodate any fMRI dataset that has been mapped to fs_LR_32k  surface space (CIFTI format, “.dtseries.nii” file with 64k vertices, 32k vertices per hemisphere when including medial wall). For more information regarding the CIFTI format, see Glasser et al. 2013 Neuroimage (https://pubmed.ncbi.nlm.nih.gov/23668970/). 

If you have issues, please email me (Chuck Lynch; cjl2007@med.cornell.edu), and I will try my best to help you. 

**Before you begin**

Check that all of the necessary dependencies are available. 

Start Matlab, open the pfm_tutorial.m script in the Matlab editor window, and run lines 5-6 to add all of the necessary dependencies to the Matlab search path.

![1](https://github.com/cjl2007/PFM-Depression/assets/46632198/41a2c422-b5bb-4b57-8a14-28c5fc6377d7)

Next, run lines 8-10 to define the number of workers for parallelization of certain procedures later on. The appropriate number for you will depend on your computing environment and resources. 

![2](https://github.com/cjl2007/PFM-Depression/assets/46632198/f5d32b11-634b-40f1-9e14-a2606f60866b)

You will also need to obtain and install Infomap (https://www.mapequation.org/infomap/#Install). You will need to define a variable called InfomapBinary during Step 4 that contains the path to the binary executable and ensure that it can be called from your Matlab environment.

Finally, you should download the example dataset from OpenNeuro.org (https://openneuro.org/datasets/ds005118/). Please note that this is a large download (~100 GB).

![image](https://github.com/cjl2007/PFM-Depression/assets/46632198/e1f55e1b-d0ad-458b-a462-24e369377a28)

**Step 1: Temporal Concatenation of fMRI data from all sessions**

The first step is to temporally concatenate all the individual denoised and fs_LR_32k surface-registered CIFTI (“.dtseries.nii”) files to obtain a single CIFTI file for the subsequent analysis.  

For context, the validity and test-retest reliability of functional connectivity measurements and their derivatives (functional network parcellations) increases rapidly with additional data per-subject. This is why PFM is often performed in individuals that have been scanned repeatedly over an extended period of time, for total scan durations of several hours or more. For example, ME01 was scanned 44 times over 10 imaging sessions (> 10 hours of functional MRI data total). In principle, however, PFM may be performed with far less data when specialized acquisitions are used to improve signal-to-noise (such as multi-echo or high field strength). 

In the pfm_tutorial.m script, run lines 14-16 to define the path to the subject’s folder and identifier.  

![3](https://github.com/cjl2007/PFM-Depression/assets/46632198/c891e457-6386-4767-8ef5-372af41c8654)

Next, run lines 23-21 to define and create the folder (PfmDir) where all the pfm related outputs will be stored. Note that for convenience, the dataset is distributed online with this folder and all the outputs generated by the subsequent steps pre-generated. Users can delete PfmDir before proceeding if they would like to start from scratch. 

![4](https://github.com/cjl2007/PFM-Depression/assets/46632198/133b69ee-b9e5-4d5c-b355-fe86d9d9b5e3)

Next, run lines 16-42 to iteratively load each of the denoised and fs_LR_32k surface-registered CIFTI (“.dtseries.nii”) files. The resting-state time courses will be extracted from the .data field of each CIFTI file, and demeaned and motion-censored (time points with framewise displacement > 0.3 mm are discarded) before being temporally concatenated into the ConcatenatedData variable. 

![5](https://github.com/cjl2007/PFM-Depression/assets/46632198/1eff35f2-f8c6-421b-8300-0def54615cb9)

Finally, run lines 51-54 to create a single CIFTI file (ConcatenatedCifti) that contains all available resting-state data for this subject.

![6](https://github.com/cjl2007/PFM-Depression/assets/46632198/262995d7-ebae-4580-a4a2-75fe1819ece7)

**Step 2: Create distance matrix & regress nearby cortical signals from subcortical structures**

The next step is to create a matrix summarizing the distance between all points in the brain. Geodesic and Euclidean space is used for cortico-cortical (vertex-to-vertex) and subcortical-cortical distance (voxel-to-vertex), respectively. 

In the pfm_tutorial.m script in the Matlab editor window, and run lines 54-59 to define the path the the subject’s fs_LR_32k midthickness (the midpoint between the white and pial surfaces) and run the pfm_make_dmat function.

![7](https://github.com/cjl2007/PFM-Depression/assets/46632198/b4fd1d72-ed4d-4f56-92b6-e0b6a0f68391)

For context, this distance matrix is used in multiple different ways during PFM. As one example, spurious coupling between subcortical voxels and adjacent cortical tissue (e.g., inflated FC between between occipital cortex and the cerebellum) can be mitigated by regressing the average time-series of cortical tissue within a specified distance from any subcortical voxel 2,3.

In the pfm_tutorial.m script in the Matlab editor window, run lines 61-66 to run the pfm_regress_adjacent_cortex function and save the resultant CIFTI file.

![8](https://github.com/cjl2007/PFM-Depression/assets/46632198/3d2b4344-5a83-4033-b53e-8823261ea4f9)

**Step 3: Apply the desired amount of spatial smoothing**

In the Matlab editor window, the user needs to specify a range of kernel sizes (in sigma) at line 77. The example range of sigma values specified below of 0.85, 1.7, and 2.55 correspond to a FWHM of 2mm, 4mm, and 6mm, respectively (FWHM ≈ 2.355 * sigma). 

![10](https://github.com/cjl2007/PFM-Depression/assets/46632198/625be66a-708c-43e4-82fe-a2240e78c2ef)

Next, run lines 79-87 to apply the specified levels of spatial smoothing to the concatenated CIFTI file with geodesic (for cortical vertices) and Euclidean (for subcortical voxels) Gaussian kernels using Connectome Workbench command line utilities. Note that this step can be considered optional, but is recommended for most datasets to improve signal-to-noise.

![9](https://github.com/cjl2007/PFM-Depression/assets/46632198/274603c6-4ac1-4453-9984-54ebc49c1174)

**Step 4: Run infomap**

The Infomap community detection algorithm (https://www.mapequation.org/infomap/) is one of the most widely used approaches for delineating functional brain networks and their boundaries in individuals. The pfm_infomap function is a wrapper that encompasses multiple steps — including creating and thresholding the functional connectivity (FC) matrix, calling the Infomap algorithm, and saving the resultant Infomap communities to a CIFTI file.

In the Matlab editor window, run lines 91-92 to load the CIFTI file that you want to run Infomap on. In the example below, we selected the CIFTI file with 2.55 sigma spatial smoothing.

![11](https://github.com/cjl2007/PFM-Depression/assets/46632198/110da348-7008-47a0-8afa-c6367126c4d5)

Next, in the Matlab editor window, the following inputs must be defined by the user at lines 94-101 in the pfm_tutorial.m script.

![12](https://github.com/cjl2007/PFM-Depression/assets/46632198/9a5cde1a-fe09-49d5-972e-d0a8e8c7cd6d)

DistanceMatrix is the path (string) to the distance matrix created earlier during Step 2. 
DistanceCutoff is the threshold (numeric, in mm) for removing short-distance correlations in the FC matrix. Correlations between nodes < DistanceCutoff from each other will be set to zero. This is done to mitigate the effects of spatial autocorrelation on the network structures identified. 
GraphDensities is a numeric vector of graph densities — the percentage of the top connections retained by each node after thresholding. So for example, in a hypothetical 1000 x 1000 FC matrix, a graph density of 0.05 (5%) means that each node will retain its top 50 strongest connections. 

By default the total number of communities identified by Infomap is data-driven, but can be controlled in part by how many connections are retained in the functional connectivity matrix after thresholding. For example, fewer communities are identified at the 5% threshold (on average, 8.28 ± 1.21) than at the 0.1% threshold (on average, 89.13 ± 8.04). We recommend running Infomap over a range of graph densities (e.g., 5% to 0.01%, as done in 4).

![infomap](https://github.com/cjl2007/PFM-Depression/assets/46632198/d69c2b2c-9e00-4cda-839f-5694ac7de378)

NumberReps is a numeric value representing the number of times the Infomap algorithm is run before selecting the best solution. 
BadVertices is an optional index of all points in the brain the user would like to omit from the analysis. For example, if the user knows that data quality is especially poor  (low tSNR or test-retest reliability) in a particular set of brain regions. Otherwise, set to [] to include all vertices and voxels.
BrainStructures is a cell array of brain structures of interest from the .brainstructurelabel field of the ConcatenatedCifti file. The FC matrix will omit nodes from brain structures not included in this variable.
InfoMapBinary is a path (string) to the Infomap binary (see https://www.mapequation.org/infomap/). Users will need to update this example path to reflect where the Infomap binary is located on their computer.

Run lines 103-104 in the pfm_tutorial.m script to run the pfm_infomap function.

![13](https://github.com/cjl2007/PFM-Depression/assets/46632198/8a65c95d-c05d-4610-9dd3-8ff0b8b4c9fd)

Optionally, users can discard pieces (contiguous patches of vertices or voxels) of communities smaller than a specified size (in mm²). This procedure acts as a spatial filter  — removing small objects without imposing additional spatial smoothing on the underlying data. The neighboring network identities will then be dilated one vertex at a time until the region is filled, as done in 5.

![14](https://github.com/cjl2007/PFM-Depression/assets/46632198/22cf2099-516c-4cff-824e-35a19498f518)

**Step 5: Algorithmic assignment of network identities to infomap communities**

The main output of pfm_infomap is a CIFTI file called “Bipartite_PhysicalCommunities.dtseries.nii”, which contains the Infomap communities obtained at each graph density (each column represents a different graph density). These community labels are arbitrary — in other words, community number 1 will not represent the same functional network in different individuals.

We will assign each Infomap community to one of 20 known functional network identities based on their spatial locations and functional connectivity. In principle, this could be accomplished manually by an expert familiar with functional network topography in individuals, but this can be time consuming and difficult to scale. To help accelerate and standardize the network identification process, we have created a semi-automated procedure for quantifying the likelihood of an Infomap community belonging to a particular functional brain network, and specifying the best match as the initial assignment. The entire procedure is implemented using the pfm_identify_networks function.

First, in the Matlab editor window, run lines 121-122 in the pfm_tutorial.m script to load the default priors for the network identification algorithm. 

![18](https://github.com/cjl2007/PFM-Depression/assets/46632198/82ce2036-f9b8-4be8-8392-e9fd61a84965)

Priors is a structure array containing the average FC (Priors.FC) and spatial locations (Priors.Spatial) of 20 functional brain networks — (Default-Parietal, Default-Anterolateral, Default-Dorsolateral, Default-Retrosplenial, Visual-Lateral, Visual-Dorsal/Ventral Stream, Visual-V1, Visual-V5, Frontoparietal, Dorsal Attention, Premotor / Dorsal Attention II, Language, Salience, Cingulo-opercular / Action-mode, Parietal memory, Auditory, Somatomotor-Hand, Somatomotor-Face, Somatomotor-Foot, Auditory, or Somato-Cognitive-Action). The code will accept other sets of priors if they are organized in the same way.

![PriorsFC](https://github.com/cjl2007/PFM-Depression/assets/46632198/06bf2607-607d-491b-bc6a-a21b374a83bb)
![PriorsSpatial](https://github.com/cjl2007/PFM-Depression/assets/46632198/2fdaf9cb-3d36-49ab-98d4-afa1a0da5c8d)

Next, in the Matlab editor window, run lines 124-127 in the pfm_tutorial.m script to define the other inputs for the pfm_identify_networks function.

![19](https://github.com/cjl2007/PFM-Depression/assets/46632198/3a268e99-96fc-4427-9f3b-fb8624b2e297)

InfomapCommunities is the CIFTI file created by  pfm_infomap. It contains the Infomap communities obtained at each graph density (each column represents a different graph density). 
Column is the column of InfomapCommunities that the user wants to assign network identities. 
Output is the name for the output file. 

![20](https://github.com/cjl2007/PFM-Depression/assets/46632198/128354b2-4b9d-40d7-ba1b-2394ee6e5d53)

**Step 6: Review algorithmic network assignments**

There are multiple outputs generated by pfm_identify_networks that we recommend reviewing carefully before proceeding. 

First, there is an .XLS sheet containing the winning and runner-up network assignments for each Infomap community, as well as a  brief summary of the information that influenced algorithmic assignments. Second, there are a series of CIFTI files that are highlighted and described below. 

![21](https://github.com/cjl2007/PFM-Depression/assets/46632198/72a22cdd-7bc7-42be-b3bd-dd8a91cecf86)

Bipartite_PhysicalCommunities+AlgorithmicLabeling.dlabel.nii is a CIFTI file containing the initial algorithmic network assignments for all Infomap communities (set in Tab 1 in images below). The colors used to represent different functional brain networks are set by the RGB values in Priors.NetworkColors.
Bipartite_PhysicalCommunities+AlgorithmicLabeling_InfoMapCommunities.dlabel.nii is a CIFTI file containing the initial algorithmic network assignment for each Infomap community, stored separately in each column (set in Tab 2 in images below). 
Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_WholeBrain.dtseries.nii is a CIFTI file containing the whole-brain functional connectivity of each Infomap community, stored separately in each column (set in Tab 3 in images below). 
Bipartite_PhysicalCommunities+AlgorithmicLabeling_FC_btwn_InfoMapCommunities.dtseries.nii is a CIFTI file containing the functional connectivity of each Infomap community to all other Infomap communities, stored separately in each column (set in Tab 4 in images below). 

To review the algorithmic network assignments, load the subject’s anatomical data and the four CIFTI files above into Connectome Workbench.  We recommend creating a 2 x 2 tile tab configuration to simultaneously view all of the relevant information.Tabs 2-4 can be set to “Yoke I” in the Overlay Toolbox, so that when you flip through the columns of one file, the others will also update. This allows you to flip through individual communities and view their algorithmic assignment and spatial locations (Tab 2), as well as their functional connectivity with the rest of the brain (Tab 3) and with other communities (Tab 4). 

For example, community 17 (highlighted below), was labeled Auditory (Column B). The spatial correlation of community 1 FC with the Auditory template FC was r = 0.88 (Column G). The average Auditory spatial probability value in community 1 was 0.57 (Column H). The total score for this winning assignment is the product of the FC_Similarity and Spatial_Score (0.88 * 0.57 = 0.50). The total score for the first runner-up assignment (Language) is also the product of the FC_Similarity and Spatial_Score (0.39 * 0.21 = 0.08). The “Confidence” of the winning assignment (0.59) is quantified in Column I, and is defined as the relative difference in total score associated with the winning and first runner-up network assignments ([0.502 -  0.083] / 0.083 = 5.07). In this case, the high confidence score and impression after visual inspection is that this is a sensible label for this community. 

![22](https://github.com/cjl2007/PFM-Depression/assets/46632198/c718360b-8588-4b38-b3ee-fa6f2f5c3d8c)
<img width="1611" alt="23" src="https://github.com/cjl2007/PFM-Depression/assets/46632198/94cb0e31-7a5a-480c-a2f3-9f0f85b0944e">

Another example of a correctly labeled community is shown below. In this case, community 19 (highlighted below), was labeled Cingulo-opercular / Action-mode.

![24](https://github.com/cjl2007/PFM-Depression/assets/46632198/8d769d92-a57e-4062-95df-b855a48f3eab)
<img width="1611" alt="25" src="https://github.com/cjl2007/PFM-Depression/assets/46632198/f60677e7-800d-425a-b615-6ca7f41540c7">

Note that this approach works well when functional networks are situated approximately in their “typical” locations. However, sometimes communities are “incorrectly” labeled, often when they have a low Spatial_Score (i.e., minimal overlap with the typical network location). These cases tend to be accompanied by low Confidence scores. Users can sort the XLS sheet by the Confidence values (as done below) and visually examine low confidence assignments. We have found that a good rule of thumb is to focus on Confidence values of 0.33 or less (highlighted below), but it is useful to carefully review all of the assignments, if time permits.

<img width="1613" alt="26" src="https://github.com/cjl2007/PFM-Depression/assets/46632198/22a8c8df-03e8-46a4-b329-d04d560c6ddc">

To change a network assignment, copy the alternative assignment name into column C (“Network_ManualDecision”). For example, community 38 (shown below)  was originally labeled Dorsal Attention but was changed to Cingulo-opercular / Action-mode. Note how Tabs 3-4 indicate strong FC with brain regions already established as belonging to Cingulo-opercular / Action-mode. 

![27](https://github.com/cjl2007/PFM-Depression/assets/46632198/bf71f0f1-a996-4ebd-97ce-d312e7a94784)
<img width="1613" alt="Screen Shot 2024-05-07 at 2 20 41 PM" src="https://github.com/cjl2007/PFM-Depression/assets/46632198/2be6c99e-362a-4b60-84a2-f67fac4c8771">


Once all of the assignments have been inspected and the XLS file containing the manual decisions has been saved, you can run pfm_parse_manual_decisions to incorporate the manual labels.

In the Matlab editor window, run lines 134-136 in the pfm_tutorial.m script to define the inputs. 

![Screen Shot 2024-05-07 at 2 22 54 PM](https://github.com/cjl2007/PFM-Depression/assets/46632198/90b2c70b-1518-440a-818e-280a0a224908)

Next, in the Matlab editor window, run lines 138-139 to run pfm_parse_manual_decisions. 

![Screen Shot 2024-05-07 at 2 23 21 PM](https://github.com/cjl2007/PFM-Depression/assets/46632198/6f3d46c6-83b8-482d-bf50-89822da044df)

**Step 7: Calculate size of each functional brain network**

The main output of pfm_parse_manual_decisions is a CIFTI file (called “Bipartite_PhysicalCommunities+FinalLabeling.dlabel.nii”). The pfm_calculate_network_size function can be used to calculate the size of each functional brain network (the percentage of cortical surface area it occupies).

In the Matlab editor window, run lines 143-146 in the pfm_tutorial.m script to define the inputs. 

![Screen Shot 2024-05-07 at 1 10 06 PM](https://github.com/cjl2007/PFM-Depression/assets/46632198/7a03123f-5880-4393-b050-6907c5fc5642)

FunctionalNetworks is a CIFTIl file containing the “final” (post manual review) functional network maps. 
VA is a CIFTI file containing the surface area (in mm^2) each vertex is responsible for. 
Structures is a cell array of brain structures. The network size calculation will be constrained to the structures in this variable.

In the Matlab editor window, run lines 148-173 to run the pfm_calculate_network_size function and visualize the results.

![Screen Shot 2024-05-07 at 1 14 42 PM](https://github.com/cjl2007/PFM-Depression/assets/46632198/407fa49a-7dcf-435a-9646-8fd66cc16830)
<img width="1121" alt="Screen Shot 2024-05-07 at 2 26 48 PM" src="https://github.com/cjl2007/PFM-Depression/assets/46632198/06e6c847-6faa-4310-97dd-37658757cd18">
