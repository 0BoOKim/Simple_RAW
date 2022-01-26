clear all;
close all;

BEB_ENABLE = 0;
    
%-------------------------------------------
% ASSUMPTIONS
%-------------------------------------------
% A1. Carrier sensing is perfect so that all the users can completely 
% detect busy channel
% A2. There is no transmssion failure due to poor quality of wireless
% channel, i.e., transmission failure only occurs due to collision
% A3. All the users always have data to send 

%------------------------------------------
% PARAMETERS & VARIABLES
%------------------------------------------
N_slot = 1000000;         % total number of slots
N_user = 6; %2/10/50/100/200            % number of users

R_data = 0.3*10^6*ones(1,N_user); 
%R_data = [12 24 48]*10^6;
                        % transmission rate 

CW_min = 16*ones(1,N_user);             % minimum contention window
%CW_min = [64 32 16];

AIFSN = 3*ones(1,N_user);  % AIFSN for EDCA, default = 3 considering DIFS
%AIFSN = [4 4 4];

L_data = 100*ones(1,N_user);
                        % data frame size (byte)
%L_data = [1000 1000 1000];

% 802.11g MAC/PHY spec.  -->    ('22.01.12) 
T_slot = 52*10^-6;            % time slot = 20 us
L_macH = 28;                 % length of the MAC Header (bytes)
T_phyH = 44*10^-6;           % PHY Header transmission time (sec)
T_sifs = 16*10^-6;           % SIFS time (sec)
L_ack  = 14;                 % length of ACK (byte)
R_ack  = 0.3*10^6;%6*10^6;    % basic rate for ACK (bit/sec)
T_ack  = T_phyH + L_ack*8/R_ack;    
                        % ACK transmission time (sec)
T_data = T_phyH + (L_macH + L_data)*8./R_data; 
                        % data frame transmission time (sec)

T_txslot = ceil( (T_data + T_sifs + T_ack)/T_slot);   
            % number of slots required to transmit one data frame

CW_max = 1024;              
            % when BEB is enabled, CW is doubled if collision occurs
            % its maximum value is CW_max
CW = CW_min.*ones(1,N_user);                % initial contention window 

% RAW setting

On_RAW = 1;
if On_RAW
    On_CSB = 1; 
    % Option for Cross Slot Boundary (CSB) 
    % 0: CSB is disabled.
    % 1: CSB is enabled.

    % L_Beacon = 0;  % length of Beacon Frame. (TBD)
    % % Beacon_Interval = 100*10^-6; % (sec) beacon interval
    % % Beacon_Interval = ceil(Beacon_Interval/T_slot);




    C = 500;       % Slot Duration Count
    T_RAW_Slot = 500*10^-6 + C*120*10^-6;  % Duration of RAW Slot 
    N_RAW = 3;   % number of RAW Slots

    T_RAW_Slot_in_Slot = ceil(T_RAW_Slot/T_slot);
    L_Beacon = 0;  % length of Beacon Frame. (TBD)
    Beacon_Interval = T_RAW_Slot_in_Slot*N_RAW; % (slots) beacon interval
    Next_BTT = T_RAW_Slot_in_Slot;  % Next_Beacon_Transmission Time(BTT)

    % initialize RAW related parameters
        % Next_RAW_Duration = 2; % Beacon_Interval;
        % STATE_RAW = 0;    
%     Tab_RAW = [T_RAW_Slot_in_Slot  T_RAW_Slot_in_Slot+T_RAW_Slot_in_Slot 0];
%     if N_RAW > 1
%         for idx_RAW = 2:N_RAW
%             Tab_RAW(idx_RAW, 1) = Tab_RAW(idx_RAW-1, 2) + 1;
%             Tab_RAW(idx_RAW, 2) = Tab_RAW(idx_RAW, 1)+T_RAW_Slot_in_Slot;
%             Tab_RAW(idx_RAW, 3)  = 0;
%         end  
%     end
    N_Assigned_STA = ceil(N_user/N_RAW);  % number of assigned STAs in each RAW duration.

    if N_Assigned_STA > N_user
        error('N_Assigned_STA cannot be larger than N_user!');
    end
end


%initial backoff counter & aifs
bc = ceil(rand(1,N_user).*CW);
aifs = AIFSN;
aifs_reset = zeros(1,N_user);

tx_state = zeros(N_slot,N_user);
% tx_state(i,j) = transmission status at time slot i for user j
STATE_BC  = 0;   % backoff state
STATE_TX  = 1;   % transmission state (without collision)
STATE_CS  = 2;   % carrier-sensing state
STATE_COL = 3;   % collision state

n_txnode = 0; % number of TX nodes


n_access = zeros(1,N_user);
    % number of channel access per user
n_collision = zeros(1,N_user);
    % number of collisions per user
n_success = zeros(1,N_user);
    % number of successful channel access without collision
    % n_access = n_collision + n_success

i=2;
if On_RAW
    i = T_RAW_Slot_in_Slot;
end
%------------------------------------------------------------------
% START SIMULATION
%------------------------------------------------------------------
while(i < N_slot-1)
    
    
    if On_RAW
        if i == Next_BTT
           
            % configure the RAW slot
           
           Tab_RAW = [i i+T_RAW_Slot_in_Slot 0];
           
           if N_RAW > 1
               for idx_RAW = 2:N_RAW
                    Tab_RAW(idx_RAW, 1) = Tab_RAW(idx_RAW-1, 2) + 1;
                    Tab_RAW(idx_RAW, 2) = Tab_RAW(idx_RAW, 1)+T_RAW_Slot_in_Slot;
                    Tab_RAW(idx_RAW, 3)  = 0;
               end  
           end
           
           % allocate STAs for each RAW slot
           All_User = 1:N_user;
           List_Assigned_STA = [];
           
           for idx_RAW = 1:N_RAW
               for idx = 1:N_Assigned_STA
                   if isempty(All_User)
                      break;  
                   end
                   Temp_Assigned_STA = All_User(ceil(rand*length(All_User)));
                   List_Assigned_STA(idx_RAW, idx) = Temp_Assigned_STA;
                   All_User(All_User == Temp_Assigned_STA) = [];
                   
               end
           end
           % update Next Beacon Transmission Time
           Next_BTT = Next_BTT + Beacon_Interval;
        end
        
        
        % check if there exist any activated RAW slot
        for idx_RAW = 1:N_RAW
            if Tab_RAW(idx_RAW, 1) == i % when the start time of a RAW slot has reached
                if Tab_RAW(idx_RAW, 3) ~= 0 % if the STATE of a RAW slot is disabled here, something is wrong.
                     error('Unexpected RAW STATE!');
                end
                Tab_RAW(idx_RAW, 3) = 1; % the STATE of a RAW slot is changed to the enabled.
                
            end
            
             if Tab_RAW(idx_RAW, 2) == i % when the end time of a RAW slot has reached
                 if Tab_RAW(idx_RAW, 3) ~= 1 % if the STATE of a RAW slot is enabled here, something is wrong.
                     error('Unexpected RAW STATE!');
                 end
                 Tab_RAW(idx_RAW, 3) = 0;  % the STATE of a RAW slot is changed to the disabled.
             end
             
        end
        
    
    end
    
    i = i+1;
    
end