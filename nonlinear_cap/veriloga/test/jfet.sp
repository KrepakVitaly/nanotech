Verilog-A version of the SPICE JFET

.hdl jfet.va
.options post=1

VCC Drain 0 3.0
VG Gate 0 0.5
VS Source 0 0.0

X1 Drain Gate Source jfet Vto=-2.0 Beta=1.1e-4 Lambda=0.01

.dc VCC 0.0 4.0 0.01 VG -2.0 0.0 0.5

.print I(VCC)

.end