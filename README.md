LibForcesField
==============

Store and write volume field representations of forces and moments Same as *writeFields yes;* option in OpenFoam by ESI-Group

Use:

```
forceField
{
    type          forcesField;

    libs          ("libforcesField.so");

    writeControl    adjustableRunTime;
    writeInterval   0.001;

    patches       (<patch list>);
    rho           rhoInf;     // Indicates incompressible
    rhoInf        1;          // Redundant for incompressible

    CofR          (0 0 0);    // Rotation around centre line of propeller
    pitchAxis     (0 1 0);patches         (<list of patch names>);


    // Optional entries

    // Field names
    p               p;
    U               U;
    rho             rho;

    // Reference pressure [Pa]
    pRef            0;

    // Include porosity effects?
    porosity        no;

    // Store and write volume field representations of forces and moments
    writeFields     yes;

    // Centre of rotation for moment calculations
    CofR            (0 0 0);
}
```


