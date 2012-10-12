%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PDQ README
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

---------------------
      Contents
---------------------
1) Version History
2) Compatibility
3) Dependencies
4) PDQ Process
5) PDQ Function List
6) PDQ Object: MRIdata
7) PDQ Object: Template




-----------------
1) Version History
-----------------

PDQ-8-27-12
About half of dependencies added, documentation still sparse.

PDQ-6-1-11
Original posting of code, no dependencies included, sparse documentation.





-----------------
2) Compatibility 
-----------------
PDQ runs in MATLAB and has been tested in:
	-Windows Vista 64-bit (no pre-processing, see note below)
	-Windows 7 64-bit (no pre-processing, see note below)
	-Mac OS X
	-Ubuntu 10.02.

Note: PDQ cannot complete MRI data pre-processing on Windows platforms. Its phase unwrapping pre-processing step depends on a binary executable that runs only in Mac OS X and Linux. 





-----------------
0) Dependencies
-----------------
PDQ is dependent on the following MATLAB functions created by 3rd parties:

(Provided)	localmaxmin.m	(D.C. Hanselman, University of Maine, Orono, ME 04469)
(Provided)	dsxy2figxy.m	(SigPlot, Yannis Agiomyrgiannakis)
(Provided)	Simplex.m 	(Fred Sigworth
)
(Not provided)	load_nii.m	(NIFTI, XXXXXXXXX)
		save_nii.m		"
(Not provided)	normxcorr3.m	(XXXXXXXX)
		normxcorr2.m	        "
(Provided)	prelude		(FSL, XXXXXXXXXXXXXX)




-------------------------------------
1) PDQ Process
-------------------------------------

A) Import data
MRIdata = import_fid(filename, k_space_type, ordering, RARE_factor, read_dim, reco_dim, noise, NEX_scheme_y,NEX_scheme_z, bigEndian);

B) Process data
MRIdata = PDQ(MRIdata, dipole_orientation, expected_radius)


If on UNIX, you're finished. If on Windows, continue

C) Unwrap data using PRELUDE in UNIX

D) Import unwrapped data
MRIdata = prelude_import(MRIdata, path);

D) Continue processing
MRIdata = PDQ(MRIdata, dipole_orientation, expected_radius)





-------------------------------------
2) PDQ Function List
-------------------------------------




-------------------------------------
2) PDQ Data Structure: MRIdata
-------------------------------------
          .fid            (3D Float)         K-space (Original MRI acquisition echos)
          .complex        (3D Float)         Complex image
          .mag            (3D Float)         Magnitude image
          .phase          (3D Float)         Phase image
          .phase_hamming  (3D Float)         Hamming-filtered phase
          .mask           (3D Logical)       Image mask

          .unwrapped           (3D Float)    Unwrapped phase image
          .unwrapped_hamming   (3D Float)    Unwrapped hamming-filtered phase
          .rauscher            (3D Float)    Rauscher phase image (unwrapped - unwrapped_hamming)
          .ramp                (3D Float)    Phase-ramp-removed unwrapped phase
          .high_pass           (3D Float)    High-pass filtered unwrapped phase


          .PDQ                 Results from a run of PDQ - Phase Map Cross-Correlation Detection and Quantification

	      .dipoles_ramp                Array of all dipoles found in ramp-removed phase image
              .dipoles_hp                  Array of all dipoles found in high-passed phase image
              .dipoles_hp.index            (Float)    Dipole's index
              .dipoles_hp.x                (Float)    Dipole's x coordinate in original phase_unwrapped_dataset
              .dipoles_hp.y                (Float)    Dipole's y coordinate in original phase_unwrapped_dataset
              .dipoles_hp.z                (Float)    Dipole's z coordinate in original phase_unwrapped_dataset
              .dipoles_hp.phase            (3D Float) Dipole's 3D phase impression extracted from original phase_unwrapped_dataset
              .dipoles_hp.xcorr            (Float)    Dipole's XCORR value against template
              .dipoles_hp.suscept_LSF      (Float)    Dipole's Least-Squares-Fit multiple of template
              .dipoles_hp.neighbors        (1D Float) Array of dipole's neighbors (dipoles within template-sized region around dipole)
              .dipoles_hp.num_neighbors    (Float)    The number of neighbors that a dipole has (i.e., length of neighbors array)
              .dipoles_hp.cluster2         (Float)    Cluster this dipole is assigned when performing 2-Cluster Mahalanobis Furthest-Distance Clustering
              .dipoles_hp.cluster3         (Float)    Cluster this dipole is assigned when performing 3-Cluster Mahalanobis Furthest-Distance Clustering
              .dipoles_hp.cluster4         (Float)    Cluster this dipole is assigned when performing 4-Cluster Mahalanobis Furthest-Distance Clustering
              .dipoles_hp.cluster5         (Float)    Cluster this dipole is assigned when performing 5-Cluster Mahalanobis Furthest-Distance Clustering

              .neighbor_consensus       (3D Float) Consensus dipoles for different numbers of neighbors
           

          .PSM                 Results from a run of PSM - Phase Slope Magnitude
              .PSM.mag                 (3D Float)    Phase slope magnitude image   (How rapidly phase is changing in space)
              .PSM.x                   (3D Float)    X-direction Phase slope magnitude image
              .PSM.y                   (3D Float)    Y-direction Phase slope magnitude image
              .PSM.z                   (3D Float)    Z-direction Phase slope magnitude image
              .PSA                     (3D Float)    Phase slope angle image       (The angle (in radians) at which phase is changing)
              .AC                      (3D Float)    Phase angle coherence image   (Measure of how regionally consistent PSA image is)
              .Composite               (3D Float)    Composite image               (PSM * AC. Highlights commonly-oriented rapid phase changes!)


          .regions             Region-of-interest data
                  .mask                (3D Logical)  Region mask
                  .dipole_count        (Float)       Number of dipoles contained in the region
                  .dipole_min_dist     (1D Float)    List of minimum distances of dipoles from this region





-------------------------------------
3) PDQ Data Structure: Template
-------------------------------------
           .template       (3D Float)    Theoretical template phase data impression
           .resolution     (1D Float)    Template pixel resolution (microns) [x_res y_res z_res]
           .a              (Float)       Template dipole radius, "a"
           .B0             (Float)       Primary magnetic field (MHz)
           .TE             (Float)       Echo time (milliseconds)
           .CHI_spheroid   (Float)       Volume magnetic susceptibility of spheroid's material (unitless)
           .CHI_media      (Float)       Volume magnetic susceptibility of background medium in which sphere is suspended (unitless)
           .model_inside   (Logical)     Flag (1 or 0) denoting whether template was generated modeling the inside of the sphere or not






-------------------------------------
4) xxxxx
-------------------------------------


