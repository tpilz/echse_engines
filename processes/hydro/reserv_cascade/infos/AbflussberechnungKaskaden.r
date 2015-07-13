
## Formeln für Abflussberechnung aus Kaskaden (aus Bollrich 2007)

## 1: Überfallformel nach Poleni; Annahme: Es gibt keinen Rückstau => vollkommener Überfall
##### Q = 2/3 * p * b * sqrt(2*g) * h^3/2
##### p = Überfallbeiwert, typische Größe zwischen 0.6 und 0.7 (Freimann 2009), hier Annahme 0.65 
##### b = Wehrbreite => eine Breite für jedes subcatchment ("Breite_Damm"-Spalte berechnet in Kaskaden_Data.txt)
##### h = Höhe der Wassersäule oberhalb der Staumauer = Überstauhöhe
##### hmax = Einstauhöhe für jedes subcatchment (berechnet in Kaskaden_Data.txt)

Poleni_Q=function(b,h,hmax)

{
  2/3 * 0.65 * b * sqrt(2*9.81) * (h-hmax)^3/2
}


## 2: Formel für Outlet am Fuß des Dammes: Q = p * A * sqrt(2*9.81*h1)
##### A = 0.10m*0.10m = 0.01m²
##### Abflussbeiwert p = 0.58, für Seitenverhältnis von A = 1:1
##### h als Höhenlage des Schwerpunktes der Austrittsfläche A
##### Formel gilt nur für h/A > 1 
     
Outlet_Q=function(h)

{  
  0.58 * 0.01 * sqrt(2*9.81*(h-0.05))
}
