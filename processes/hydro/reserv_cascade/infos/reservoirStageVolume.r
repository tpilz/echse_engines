rm(list=ls())

# Benetzte Querprofilfläche als Funktion des Wasserstands, a(h)
afh= function(h,kf,ke) { kf*(h^ke) }

# Volumen als Funktion des Wasserstands, v(h)
# > Erhalten durch Integration von a(h) über Speicherlänge bei geg. Gefälle
vfh= function(h,kf,ke,s) { kf/(ke+1)/s*h^(ke+1) }

# Wasserstand als Funktion des Volumens, h(v)
# > Erhalten als Umkehrfunktion von v(h)
hfv= function(v,kf,ke,s) { (v*(ke+1)*s/kf)^(1/(ke+1)) }

# Test-Parameter
kf= 1.78    # Parameter in Funktion a(h)=kf*h^ke
ke= 2.09    # Parameter in Funktion a(h)=kf*h^ke
s= 1/10     # Gefälle in m/m

# Ausgabe von h(v) für Wasserstände an Sperrstelle bis zu einem Bsp.-Wert hmax
hmax= 5
vmax= vfh(hmax,kf,ke,s)
plot(c(0,vmax), c(0,hmax), type="n", xlab="v", ylab="h")
for (v in seq(0, vmax, length.out=100)) {
  points(v, hfv(v,kf,ke,s))
}

# Test der Lösung durch ganz billige numerische Methode; nur für hmax
lmax= hmax/s
n= 1000
dx= lmax/n
v= 0
for (i in 1:n) {
  v= v + dx * afh(s*(i*dx-dx/2),kf,ke)
}
print(paste("v bei hmax (numerisch):",v))
print(paste("v bei hmax (analytisch):",vmax))

