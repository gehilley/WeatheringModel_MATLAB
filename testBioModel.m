% % testBioModel.m
% 
% BM = bioModel(100,0.1,0.3,0);
% BM.addSolidComponent('Na',1,0);
% %BM.addSolidComponent('Ca',5,0);
% BM.addFluidComponent('NaAq',0,0);
% %BM.addFluidComponent('CaAq',0,1);
% 
% BM.addTransferFunction('Na->NaAq','Na','NaAq',1,1,1e-2);
% BM.addTransferFunction('NaAq->Na','NaAq','Na',1,1,1);
% 
% %BM.addTransferFunction('Ca->none','Ca','NULL',1,1,1e-2);
% 
% %BM.addScalingFunctionToTransferFunction('Na->NaAq','Ca',1,2,3,4,5,6);
% 
% BM1 = BM;
% 
% for(i=1:150*50)
%     
%    BM1 = BM1.integrateToTime(1);
% 
%    BM1.plotComponents
%    subplot(1,2,1)
%    hold on
%    subplot(1,2,2)
%    hold on
%    drawnow
%    
% end

BM = bioModel(100,0.1,0.3,0);
BM.addSolidComponent('FeII',1,0);
%BM.addSolidComponent('Ca',5,0);
BM.addFluidComponent('NaAq',0,0);
%BM.addFluidComponent('CaAq',0,1);

BM.addTransferFunction('Na->NaAq','Na','NaAq',1,1,1e-2);
BM.addTransferFunction('NaAq->Na','NaAq','Na',1,1,1);

%BM.addTransferFunction('Ca->none','Ca','NULL',1,1,1e-2);

%BM.addScalingFunctionToTransferFunction('Na->NaAq','Ca',1,2,3,4,5,6);

BM1 = BM;

for(i=1:150*50)
    
   BM1 = BM1.integrateToTime(1);

   BM1.plotComponents
   subplot(1,2,1)
   hold on
   subplot(1,2,2)
   hold on
   drawnow
   
end

