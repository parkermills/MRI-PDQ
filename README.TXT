%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PDQ README
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





********* DRAFT ** DRAFT ** DRAFT *************




-----------------------
      Contents
-----------------------
0) Software Dependencies
1) PDQ Process
2) PDQ Function List
2) PDQ Data Structure: MRIdata
3) PDQ Data Structure: Template
4) xxxxx



-----------------------
0) Software Dependencies
-----------------------
PDQ is dependent on MATLAB functions created by 3rd parties.
1) localmaxmin.m
2) dsxy2figxy.m
3) Simplex.m
4) NIFTI (load_nii.m / save_nii.m / etc.)
5) normxcorr3.m / normxcorr2.m




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


