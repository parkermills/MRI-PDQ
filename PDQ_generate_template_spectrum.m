%%%%%%%%%%%%%%%%%%%%%%%
%
% Template spectrum is created using radii that span from 0 microns to (max(resolution)) * sqrt(6) microns
% For example, if the maximum voxel size is 150 microns, then the maxiumum radius in the template_spectrum will be 367.4 microns.
% Sqrt(6) is used as an irrational number so that no aliasing artifacts appear in the template-generation process.
% In teneral, this maximum value is to ensure that no dipole is so large that it grows outside the size of the template
% 
%%%%%%%%%%%%%%%%%%%%%%%

function template_spectrum = PDQ_generate_template_spectrum(sample_template, num_of_radii, radius_range, template_shift)

%% Initialize variables
template_spectrum = struct(sample_template);



%% If orientation not symmetrical about Z-axis, compute all dipole shifts
if(sample_template.orientation ~= 2)
    for z_shift = 1:3
        for y_shift = 1:3
            for x_shift = 1:3
                parfor i = 1:num_of_radii
                    template_spectrum(x_shift,y_shift,z_shift,i) = PDQ_generate_template(sample_template.resolution, sample_template.orientation, radius_range(1) + i * (radius_range(2)-radius_range(1)) / num_of_radii, sample_template.B0, sample_template.TE, 1.0e-5, [template_shift*(x_shift-2) template_shift*(y_shift-2)  template_shift*(z_shift-2)]);
                end
            end
        end
    end
    
%% If orientation is symmetrical about Z-axis, take advantage and duplicate already-generated templates
else
    for z_shift = 1:3
        for y_shift = 1:2
            for x_shift = y_shift:2
                parfor i = 1:num_of_radii
                    template_spectrum(x_shift,y_shift,z_shift,i) = PDQ_generate_template(sample_template.resolution, sample_template.orientation, radius_range(1) + i * (radius_range(2)-radius_range(1)) / num_of_radii, sample_template.B0, sample_template.TE, 1.0e-5, [template_shift*(x_shift-2) template_shift*(y_shift-2)  template_shift*(z_shift-2)]);
                end
            end
        end
    end
    
    for i = 1:num_of_radii
        % Copy symmetrical templates (there exist two others)
        %template_spectrum(1,1,1,i) %top corner
        template_spectrum(1,3,1,i) = template_spectrum(1,1,1,i);
        template_spectrum(3,3,1,i) = template_spectrum(1,1,1,i);
        template_spectrum(3,1,1,i) = template_spectrum(1,1,1,i);
        for j = 1:7
            template_spectrum(1,3,1,i).template(:,:,j) = rot90(template_spectrum(1,1,1,i).template(:,:,j));
            template_spectrum(3,3,1,i).template(:,:,j) = rot90(rot90(template_spectrum(1,1,1,i).template(:,:,j)));
            template_spectrum(3,1,1,i).template(:,:,j) = rot90(rot90(rot90(template_spectrum(1,1,1,i).template(:,:,j))));
        end
        
        %template_spectrum(2,1,1,i) %top side
        template_spectrum(1,2,1,i) = template_spectrum(2,1,1,i);
        template_spectrum(3,2,1,i) = template_spectrum(2,1,1,i);
        template_spectrum(2,3,1,i) = template_spectrum(2,1,1,i);
        for j = 1:7
            template_spectrum(1,2,1,i).template(:,:,j) = rot90(template_spectrum(2,1,1,i).template(:,:,j));
            template_spectrum(3,2,1,i).template(:,:,j) = rot90(rot90(template_spectrum(2,1,1,i).template(:,:,j)));
            template_spectrum(2,3,1,i).template(:,:,j) = rot90(rot90(rot90(template_spectrum(2,1,1,i).template(:,:,j))));
        end
        
        %template_spectrum(1,1,2,i) %center corner
        template_spectrum(1,3,2,i) = template_spectrum(1,1,2,i);
        template_spectrum(3,3,2,i) = template_spectrum(1,1,2,i);
        template_spectrum(3,1,2,i) = template_spectrum(1,1,2,i);
        for j = 1:7
            template_spectrum(1,3,2,i).template(:,:,j) = rot90(template_spectrum(1,1,2,i).template(:,:,j));
            template_spectrum(3,3,2,i).template(:,:,j) = rot90(rot90(template_spectrum(1,1,2,i).template(:,:,j)));
            template_spectrum(3,1,2,i).template(:,:,j) = rot90(rot90(rot90(template_spectrum(1,1,2,i).template(:,:,j))));
        end
        
        %template_spectrum(2,1,2,i) %center side
        template_spectrum(1,2,2,i) = template_spectrum(2,1,2,i);
        template_spectrum(3,2,2,i) = template_spectrum(2,1,2,i);
        template_spectrum(2,3,2,i) = template_spectrum(2,1,2,i);
        for j = 1:7
            template_spectrum(1,2,2,i).template(:,:,j) = rot90(template_spectrum(2,1,2,i).template(:,:,j));
            template_spectrum(3,2,2,i).template(:,:,j) = rot90(rot90(template_spectrum(2,1,2,i).template(:,:,j)));
            template_spectrum(2,3,2,i).template(:,:,j) = rot90(rot90(rot90(template_spectrum(2,1,2,i).template(:,:,j))));
        end
        
        %template_spectrum(1,1,3,i) %bottom corner
        template_spectrum(1,3,3,i) = template_spectrum(1,1,3,i);
        template_spectrum(3,3,3,i) = template_spectrum(1,1,3,i);
        template_spectrum(3,1,3,i) = template_spectrum(1,1,3,i);
        for j = 1:7
            template_spectrum(1,3,3,i).template(:,:,j) = rot90(template_spectrum(1,1,3,i).template(:,:,j));
            template_spectrum(3,3,3,i).template(:,:,j) = rot90(rot90(template_spectrum(1,1,3,i).template(:,:,j)));
            template_spectrum(3,1,3,i).template(:,:,j) = rot90(rot90(rot90(template_spectrum(1,1,3,i).template(:,:,j))));
        end
        
        %template_spectrum(2,1,3,i) %bottom side
        template_spectrum(1,2,3,i) = template_spectrum(2,1,3,i);
        template_spectrum(3,2,3,i) = template_spectrum(2,1,3,i);
        template_spectrum(2,3,3,i) = template_spectrum(2,1,3,i);
        for j = 1:7
            template_spectrum(1,2,3,i).template(:,:,j) = rot90(template_spectrum(2,1,3,i).template(:,:,j));
            template_spectrum(3,2,3,i).template(:,:,j) = rot90(rot90(template_spectrum(2,1,3,i).template(:,:,j)));
            template_spectrum(2,3,3,i).template(:,:,j) = rot90(rot90(rot90(template_spectrum(2,1,3,i).template(:,:,j))));
        end
    end
end

save template_spectrum template_spectrum -v6 % In case PDQ process fails, make sure to export a copy of the template_spectrum

%%%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%%%