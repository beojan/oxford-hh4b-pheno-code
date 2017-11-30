-- Writes a test Ntuple file according to the standard format
local outfile = "test_data.dat"
io.output(outfile)
io.write("# signal source weight test_variable\n")
local n_events = 1000
for i=1, n_events, 1 do
	io.write('1 SIGNAL ', torch.uniform(),' ', torch.normal(1), '\n') -- Signal events have 'test kinematic' as a gaussian around 1
	io.write('0 BACKGR ', torch.uniform(),' ', torch.normal(0), '\n') -- Background events have 'test kinematic' as a gaussian around 0
end
