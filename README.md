libMomentsField
==============
Calculates the moments by integrating the pressure and skin-friction forces over a given list of patches.

Store and write volume field representations of moments. Same as *writeFields yes;* option in [OpenFoam by ESI-Group](https://www.openfoam.com/documentation/cpp-guide/html/guide-fos-forces-forces.html) but only to write the fields. If you want the force log file you can combine it with the classic libforce function.

Example of function object specification:
```
forceField_1
{
    type          forcesField;
    libs          ("libforcesField.so");
    ...
    patches       (<patch list>);
}
```

|Property     | Description             | Required    | Default value|
| :---------- | :---------------------- | :---------: | :----------: |
|type         | Type name: forcesField  | yes         |              |
|patches      | Patches included in the forces calculation | yes |   |
|p            | Pressure field name     | no          | p            |
|U            | Velocity field name     | no          | U            |
|rho          | Density field name (see below) | yes   | rho          |
|CofR         | Centre of rotation (see below) | no   |              |
|directForceDensity | Force density supplied directly (see below)|no|no|
|fD           | Name of force density field (see below) | no | fD     |

Note:
- For incompressible cases, set `rho` to `rhoInf`.  You will then be required to provide a `rhoInf` value corresponding to the free-stream constant density.
- If the force density is supplied directly, set the `directForceDensity` flag to `yes`, and supply the force density field using the `fDName` entry.
- The centre of rotation (`CofR`) for moment calculations can either be specified by an `CofR` entry, or be taken from origin of the local coordinate system. For example:

```
    CofR        (0 0 0);
```

or

```
    origin  (0 0 0);
    coordinateRotation
    {
        type    axesRotation;
        e1  (1 0 0);
        e2  (0 1 0);
    } 
```
