![img](pyTANSPEC_logo.png)

[Home](Pipeline_Documentation.html)

[Pipeline In a Nutshell](Pipeline_in_a_nutshell.html)

[Installation](Installation.html)

[Getting Started](Getting_started.html)

[Tasks](Tasks.html)

[Advantages and Limitations of the pipeline](adv_disadv.html)

[Appendix](Appendix.html)


# Advantages and Limitations of the pipeline


## Advantages

-   The pipeline can reduce both LR and XD spectra.
-   The pipeline can do wavelength canibration for all the slits.
-   Each task is independent to each other. So user can do any task themselves and avoid that task to be done by the pipeline.
-   Pipeline will save each intermediate file. This make it easy to debug any issues.
-   The pipeline is running in terminal. So the user can learn and handle easily.
-   The pipeline uses template matching for wavelength calibration. So the wavelength calibration will be more faster and precise than the previous version.
-   The pipeline can do flux calibration for both LR and XD mode.
-   This pipeline can be used with other NIR and Optical instruments also, by adding the trace and template files.


## Disadvantages

-   Telluric corrections are not incorporated.
-   This version does not provide the optimal extraction.

