
for i=1:20
Frame = FRAME_NODE{i};
if i < 10
    fn = ['E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\Database\UniformFlood_MilanLongNokill_Nd0', int2str(i),'.mat'];
else
    fn = ['E:\TestLRF\PERA_PlanarReconstructionAlgorithm\PeraScripts\Database\UniformFlood_MilanLongNokill_Nd', int2str(i),'.mat'];

end
save(fn, 'Frame');
end