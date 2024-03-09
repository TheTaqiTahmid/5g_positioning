clear
close all;
clc
warning('off');
% Physical constant
c = physconst('LightSpeed');

% Physical layer parameters
fc = 30e9;                      % Carrier frequency [Hz]
f0 = 120e3;                     % Subcarrier spacing [Hz]
n_subcarriers = 158;            % Number of subcarriers  % 1638, 1584     
BW = f0*n_subcarriers;          % Active bandwidth [Hz]

% Transmit power and antenna gains [dBm]
P_tx_plus_gains = 10;  

% Transmission direction used in PEB (downlink or uplink). However, if the
% downlink is selected, then the orientation of the MT plays a great role
% and needs to be adjusted. In uplink, BS orientations are important.
transmission_direction = 'uplink';        

% Scenario (choose one from the set {'ToA','AoA','AoA+ToA'})
scenario = ('AoA');

% ---------------------------------------------------------------------------------------
% 
% Geometry related options 
% ---------------------------------------------------------------------------------------

% Limits of the square area [0,area_limit]x[0,area_limit] 
area_limit = 400;               

% BS locations (matrix with dimensions [2 x number_of_BSs])
BS_location2 = [200,200,200,200,250,250,250,250,300,300,300,300;200,200,200,200,300,300,300,300,200,200,200,200];
BS_location = [200, 250, 300; 200, 300,200];
% BS_location = [50,200,350;350,50,350];
% This is for the example BS orientations (not used in this example)
angles_to_center = atan2(area_limit/2-BS_location(2,:),area_limit/2-BS_location(1,:));

% BS orientations roughly towards the MT (affects CRB, does not play a role in downlink)
%BS_orientation = angles_to_center + 0.5*randn(1,size(BS_location,2));
BS_orientation = [0,pi/2,pi,3*pi/2, 0,pi/2,pi,3*pi/2, 0,pi/2,pi,3*pi/2];
% BS_orientation = [2*pi/3,4*pi/3,0];
% MT orientation (affects CRB, does not play a role in uplink)
MT_orientation = pi/4;

% ---------------------------------------------------------------------------------------
% PEBs for each MT location on the map
% ---------------------------------------------------------------------------------------

% xy-coordinates for the MT locations on the map (1m grid)
[MT_X,MT_Y] = meshgrid(linspace(0,area_limit,area_limit));

% Initialize matrix for the PEB values
PEB = nan(size(MT_X));
SNR_matrix = nan(size(MT_X));
% Loop over all the MT locations
for i = 1:length(MT_X(:))
    
    % Current MT location
    MT_location = [MT_X(i); MT_Y(i)];
    
    % Distance to the BSs
    BS_MT_distances = sqrt(sum((BS_location2-MT_location).^2,1));
    
    % ---------------------------------------------------------------------------------------
    % Cramer-Rao bounds for the ToA and AoA measurements (NEEDS TO BE IMPLEMENTED!)
    % ---------------------------------------------------------------------------------------
    
    % SNR that takes into account the free-space path-loss and thermal noise
    % over the bandwidth
    M = (n_subcarriers-1)/2;
    M = round(M);
    Mt = (M*(M+1)*(2*M+1))/2;
    SNR_dB = P_tx_plus_gains - 20*log10(4*pi*BS_MT_distances/c*fc) + 174 - 10*log10(BW);
    SNR = 10.^(SNR_dB/10);
    toa_crbs = (1 ./ (8*pi^2*f0^2*SNR*Mt));
    %toa_crbs = 3e8.*sqrt(ToA_VAR); 
    L = 16;
    
%     y = abs(BS_location(2,:) - MT_Y(i));
%     x = abs(BS_location(1,:) - MT_X(i));
%     phi(1) = atan(y(1)/x(1))+ BS_orientation(1);
%     phi(2) = atan(y(2)/x(2))+ BS_orientation(2);
%     phi(3) = atan(y(3)/x(3))+ BS_orientation(3);
    phi_global = atan2(BS_location2(2,:)-MT_Y(i),BS_location2(1,:)-MT_X(i)); 
    phi = phi_global-BS_orientation;
    
    for jj = 1:3
        if phi(1) >= phi(2) && phi(1)>= phi(3)
            phi_act(1) = phi(1);
        elseif phi(2) >= phi(1) && phi(2) >= phi(3)
            phi_act(1) = phi(2);
        elseif phi(3) > phi(2) && phi(3)>phi(1)
            phi_act(1) = phi(3);
        end
    end
    for kk = 1:3
        if phi(4) >= phi(5) && phi(4)>=phi(6)
            phi_act(2) = phi(4);
        elseif phi(5) >= phi(4) && phi(5)>= phi(6)
            phi_act(2) = phi(5);
        elseif phi(6) > phi(5) && phi(6)>phi(4)
            phi_act(2) = phi(6);
        end
    end
    for ll = 1:3
        if phi(7) >= phi(8) && phi(7)>= phi(9)
            phi_act(3) = phi(7);
        elseif phi(8) >= phi(7) && phi(8) >= phi(9)
            phi_act(3) = phi(8);
        elseif phi(9) > phi(7) && phi(9)>phi(8)
            phi_act(3) = phi(9);
        end
    end
    
    
    kd = pi;  
    aoa_crbs = 6./(L*(L^2-1)*SNR(:,1:3:end)*kd^2.*(cos(phi_act)).^2);
    %aoa_crbs = (180/pi).*sqrt(AoA_VAR);
     
    % ToA and AoA (values below are from the "Positioning in Wireless
    % Communications Systems" for Fig. 4.7)
    %toa_crbs = ones(1,length(BS_location))/c^2;
    %aoa_crbs = 0.01*ones(1,length(BS_location))/180*pi;
    
    % ---------------------------------------------------------------------------------------
    % Position errors bounds based on the geometry and obtained AoA and ToA CRBs
    % ---------------------------------------------------------------------------------------
    SNR_matrix(i) = max(SNR_dB);
    PEB(i) = get_PEB(MT_location,BS_location,toa_crbs,aoa_crbs,scenario,transmission_direction);
    
end


% ---------------------------------------------------------


%------------------------------
% Example plotting 
% ---------------------------------------------------------------------------------------
figure;
histogram(SNR_matrix);
% In order to match the plot with Fig. 4.7, I removed the PEB values > 3m
PEB(PEB>3) = nan;

% Rotation matrix for plotting the BS orientations
R = @(a) [cos(a), -sin(a); sin(a), cos(a)];

figure;
contourf(MT_X,MT_Y,reshape(PEB,size(MT_X)),20,'EdgeColor','none');
hold on;
h_bs = plot(BS_location(1,:),BS_location(2,:),'sr','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',10);
for bs = 1:size(BS_location2,2)
    orientation_vector = R(BS_orientation(bs))*[area_limit/25;0];
    plot(BS_location2(1,bs)+[0,orientation_vector(1)],BS_location2(2,bs)+[0,orientation_vector(2)],'-k','LineWidth',2);
end
hold off;
cbar = colorbar;
set(gca,'YLim',[0,area_limit],'XLim',[0,area_limit],'FontSize',14);
set(cbar,'FontSize',14);
set(cbar.Label,'String','PEB [m]');
ylabel('y [m]');
xlabel('x [m]');
%figLabel = {'Carrier Frequency:3.5 GHz'};
%dim = [0, 0.07, 0, 0];
%annotation('textbox', dim, 'String', figLabel, 'FitBoxToText', 'on', 'LineStyle', 'none');
title('Effect of secotored BSs on AOA, Cell per BS: 3')
grid on; box on;
legend(h_bs,'BS locations','FontSize',14);
