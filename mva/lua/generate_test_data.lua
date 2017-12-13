-- Writes a test Ntuple file according to the standard format
local outfile = "test_data.dat"
io.output(outfile)
io.write("# signal source weight test_variable_1 test_variable_2\n")
local n_events = 1000
for i=1, n_events, 1 do
	io.write('1 SIGNAL ', torch.uniform(),' ', torch.normal(1), ' ', torch.normal(0), '\n') 
	io.write('0 BACKGR ', torch.uniform(),' ', torch.normal(0), ' ', torch.normal(1), '\n')
end
