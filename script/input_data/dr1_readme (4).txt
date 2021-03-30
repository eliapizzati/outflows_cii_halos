###Catalogs###

Catalogs are there. Readme in the directory.

cii_target_1_arcsec.csv = flux of the detected [CII] target sources

target_cii_upper_limits.txt = upper limits of the non detections (for targets only)

continuum_target_1_arcsec.csv = flux of the continuum-detected target sources

target_continuum_upper_limits.txt = upper limits of the non detections (for targets only)

continuum_non_target.csv = catalog of non-target continuum detections

continuum_non_target_decontaminated.csv = same as previous but after re-imaging the continuum of the non-target sources contaminated by a line using only the clean continuum channels.                  

###cii_mom0_maps (and cii_mom0_maps_tapered)###

Contain the moment-0 maps of [CII].  The channels used for this moment map are the one where the target source is emitting (bright non-target sources at other frequencies might not be in this maps). File names are explicit. The image, flux, and psf files contain the science images, the primary gain maps, and the synthesized beam maps, resp. The _tapered directory contains similar data, but with a 1.5-arcsec tapering.

### continuum_maps, continuum_maps_tapered, continuum_maps_non_target_decontaminated ###

Contain the continuum maps (tapered and not). Naming conventions are similar as for mom-maps. The channels contaminated by the [CII] line of the targets were excluded.

continuum_maps_non_target_decontaminated contains special maps in which we re-measured the continuum of non-target sources after excluding the channels contaminated by a bright line (CO or [CII]).

### raw_cubes ###

Contains the raw cubes (continuum is included)
