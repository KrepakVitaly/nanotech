Verilog-A My_Model

.hdl my_model.va
*.options post=0

.probe V(port1,0)
.probe V(0)
*.print V(port1)

VCC Drain 0 3.0
VG Gate 0 0.5
VS Source 0 0.0

*X1 Drain Gate Source jfet Vto=-2.0 Beta=1.1e-4 Lambda=0.01

X1 port1 0 my_model

.tran 100n 10m

*.dc VCC 0.0 4.0 0.01 VG -2.0 0.0 0.5

*.print I(VCC)

.end