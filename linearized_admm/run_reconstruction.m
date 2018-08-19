function run_reconstruction(scene)
% Function takes as input an integer between 1 and 5, 
% corresponding to the following scenes:
%   1 - mannequin
%   2 - exit sign
%   3 - "SU" scene 
%   4 - outdoor "S"
%   5 - diffuse "S"

scene_args = getSceneArgs(scene);     
ladmm(scene_args);       

end

        


