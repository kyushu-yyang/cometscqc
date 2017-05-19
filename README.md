# Quench Code with Propagator

---

# Author

* Ye YANG (Kyushu University)
* Email: kanouyou@kune2a.nucl.kyushu-u.ac.jp

# Overview
QuenchCode is a Toolkit for calculation of magnet quench with quench propagator.
The code is wriiten in C++ and depends on the libraries:

- ROOT
- matplotlib


# Code Sample for calculation of the magnetic field
Magnetic field of coil is calculated by using the Biot-Savart Law for solenoids.
The simple way to obtain the magnetic field for coil is to use the `XFieldHandle` class.

```cpp
std::vector<Quench::XFieldContainer*> cs1;
XFieldHandle* fld = new XFieldHandle();

// construct magnet
fld->AddCoil("CS0", 857.88*mm, 1038.12*mm, 672.*mm, 823.65*mm);
fld->SetMesh("CS0", 35, 18);

fld->AddCoil("CS1", -595.25*mm, 795.25*mm, 672.*mm, 823.65*mm);
fld->SetMesh("CS1", 270, 18);

fld->AddCoil("MS1", -2121.375*mm, -653.625*mm, 672.*mm, 756.25*mm);
fld->SetMesh("MS1", 285, 10);

fld->SetTarget("CS1");
fld->Run();
cs1 = fld->GetFieldCollection();
```

# Demo for Coil Construction
```cpp
void SetConductor(XCoilHandle* hand)
{
  XCoilConductor* cdt = new XCoilConductor();
  cdt->SetDimension(4.73*mm, 15.*mm);
  cdt->SetInsSize(0.3*mm, 0.3*mm);
  hand->SetConductor( dynamic_cast<XCoilBase*>(cdt) );
}

void SetStrip(XCoilHandle* hand)
{
  XCoilStrip* strip = new XCoilStrip();
  strip->SetDimension(4.73*mm+0.3*2*mm, 1.*mm);
  strip->SetInsSize(0.0, 0.3*mm);
  hand->SetStrip( dynamic_cast<XCoilBase*>(strip) );
}

void Construct()
{
  XCoilHandle* hand = new XCoilHandle();
  try {
    hand->SetName("CS1");
    hand->SetCoilSize(0., 0.672*m, 0.);
    hand->SetMesh(270, 4, 19);
    hand->SetCoilLayers(9);
    hand->SetCoilTurns(270);
    hand->SetMaterialRatio(7.3, 1., 0.9);
    SetConductor(hand);
    SetStrip(hand);
    hand->AddLayer(1, kStrip);
    hand->AddLayer(2, kConductor);
  }
  catch (XQuenchExcept except) {
    delete hand;
  }
}
```
