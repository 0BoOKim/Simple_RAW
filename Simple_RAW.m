%--------------------------------------------------------------------------
% Simple RAW (Restricted Access Window) in IEEE 802.11ah
%--------------------------------------------------------------------------
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
           
           %initial backoff counter & aifs
            bc = ceil(rand(1,N_user).*CW);
            aifs = AIFSN;
            aifs_reset = zeros(1,N_user);
            
            
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
        
        if length(find(Tab_RAW(:,3) == 1)) > 1
            error('only one RAW slot can be enabled!');
        end
        
        Idx_Enabled_RAW = find(Tab_RAW(:,3) == 1); 
        Assigned_STA = List_Assigned_STA(Idx_Enabled_RAW,:);
        
        if ismember(0, Assigned_STA)    %
            Assigned_STA(Assigned_STA == 0) = [];
        end
        
        N_user_in_RAW_Slot = length(Assigned_STA);
        
      
        % check if channel is idle
        if (n_txnode == 0)
            for j=1:N_user_in_RAW_Slot
                if (aifs(Assigned_STA(j)) > 0)
                    aifs(Assigned_STA(j)) = aifs(Assigned_STA(j)) - 1;
                    % wait for AIFS
                end
            end
            for j=1:N_user_in_RAW_Slot
                if (aifs(Assigned_STA(j)) == 0)
                    %if (aifs_reset(j) == 1)
                        bc(Assigned_STA(j)) = bc(Assigned_STA(j)) -1; 
                    %elseif (aifs_reset(j) == 0)
                    %    bc(j) = 0;
                    %end
                end

            end
            % decrement backoff counter by 1 for users after waiting for AIFSN
        end
        % if channel is busy, do not change aifs & backoff counter 

        % check whether BC = 0
        for j=1:N_user_in_RAW_Slot
            if (bc(Assigned_STA(j)) == 0)
                tx_state(i:(i+T_txslot(Assigned_STA(j))-1),Assigned_STA(j)) = STATE_TX;
                % set sate from i to i+T_txslot-1 = STATE_TX
                n_txnode = n_txnode + 1;          
                n_access(Assigned_STA(j)) = n_access(Assigned_STA(j))+1;

                bc(Assigned_STA(j)) = ceil(rand*CW(Assigned_STA(j)));
                aifs(Assigned_STA(j)) = AIFSN(Assigned_STA(j));            
                % re-select a new random backoff & aifs
            end
        end

        % update state     
        % if channel is busy
        if (n_txnode ~= 0 )
            % if at least one node is in transmision state
            max_duration = max( T_txslot.*(tx_state(i,:)==STATE_TX) );
            for (j=1:N_user_in_RAW_Slot)
                if (tx_state(i,Assigned_STA(j)) ~= STATE_TX)
                    tx_duration = i+max_duration-1;
                    tx_state(i:tx_duration,Assigned_STA(j)) = STATE_CS;
                    % set state = carrier sensing
                    aifs(Assigned_STA(j)) = AIFSN(Assigned_STA(j));
                    %aifs_reset(j) = 1;
                    % reset aifs after sensing busy channel
                end
            end

            % check collision
            if (n_txnode == 1)
                % only one node accesses channel, i.e., no collision
                for (j=1:N_user_in_RAW_Slot)
                    if (tx_state(i,Assigned_STA(j)) == STATE_TX)
                        n_success(Assigned_STA(j)) = n_success(Assigned_STA(j)) + 1;
                        %aifs_reset(j) = 0;
                        % BEB here...
                        if (BEB_ENABLE == 1)
                            CW(Assigned_STA(j)) = CW_min(j);
                        end
                    end
                end
            elseif (n_txnode > 1)
                % more than two nodes access channel => collision
                for (j=1:N_user_in_RAW_Slot)
                    if (tx_state(i,Assigned_STA(j)) == STATE_TX)
                        tx_duration = i+T_txslot(Assigned_STA(j))-1;
                        tx_state(i:tx_duration,Assigned_STA(j)) = STATE_COL;
                        % set state = collision
                        n_collision(Assigned_STA(j)) = n_collision(Assigned_STA(j))+1;
                        %aifs_reset(j) = 1;
                        % BEB here.....
                        if (BEB_ENABLE == 1)
                            CW(Assigned_STA(j)) = min(CW(Assigned_STA(j)) * 2, CW_max); 
                        end
                    end
                end        
            end % end for collision-check

            i = i + max_duration+1;   % increase time index by T_txslot
            n_txnode = 0;
        else
            i=i+1;  % increase time index by 1 (if n_txnode = 0)        
        end
        
    else % No RAW, normal DCF 
%         % check if channel is idle
%         if (n_txnode == 0)
%             for j=1:N_user
%                 if (aifs(j) > 0)
%                     aifs(j) = aifs(j) - 1;
%                     % wait for AIFS
%                 end
%             end
%             for j=1:N_user
%                 if (aifs(j) == 0)
%                     %if (aifs_reset(j) == 1)
%                         bc(j) = bc(j) -1; 
%                     %elseif (aifs_reset(j) == 0)
%                     %    bc(j) = 0;
%                     %end
%                 end
% 
%             end
%             % decrement backoff counter by 1 for users after waiting for AIFSN
%         end
%         % if channel is busy, do not change aifs & backoff counter 
% 
%         % check whether BC = 0
%         for j=1:N_user
%             if (bc(j) == 0)
%                 tx_state(i:(i+T_txslot(j)-1),j) = STATE_TX;
%                 % set sate from i to i+T_txslot-1 = STATE_TX
%                 n_txnode = n_txnode + 1;          
%                 n_access(j) = n_access(j)+1;
% 
%                 bc(j) = ceil(rand*CW(j));
%                 aifs(j) = AIFSN(j);            
%                 % re-select a new random backoff & aifs
%             end
%         end
% 
%         % update state     
%         % if channel is busy
%         if (n_txnode ~= 0 )
%             % if at least one node is in transmision state
%             max_duration = max( T_txslot.*(tx_state(i,:)==STATE_TX) );
%             for (j=1:N_user)
%                 if (tx_state(i,j) ~= STATE_TX)
%                     tx_duration = i+max_duration-1;
%                     tx_state(i:tx_duration,j) = STATE_CS;
%                     % set state = carrier sensing
%                     aifs(j) = AIFSN(j);
%                     %aifs_reset(j) = 1;
%                     % reset aifs after sensing busy channel
%                 end
%             end
% 
%             % check collision
%             if (n_txnode == 1)
%                 % only one node accesses channel, i.e., no collision
%                 for (j=1:N_user)
%                     if (tx_state(i,j) == STATE_TX)
%                         n_success(j) = n_success(j) + 1;
%                         %aifs_reset(j) = 0;
%                         % BEB here...
%                         if (BEB_ENABLE == 1)
%                             CW(j) = max(CW(j)-32, CW_min(j));%CW_min(j);
%                         end
%                     end
%                 end
%             elseif (n_txnode > 1)
%                 % more than two nodes access channel => collision
%                 for (j=1:N_user)
%                     if (tx_state(i,j) == STATE_TX)
%                         tx_duration = i+T_txslot(j)-1;
%                         tx_state(i:tx_duration,j) = STATE_COL;
%                         % set state = collision
%                         n_collision(j) = n_collision(j)+1;
%                         %aifs_reset(j) = 1;
%                         % BEB here.....
%                         if (BEB_ENABLE == 1)
%                             CW(j) = min(CW(j) * 3, CW_max); 
%                         end
%                     end
%                 end        
%             end % end for collision-check
% 
%             i = i + max_duration+1;   % increase time index by T_txslot
%             n_txnode = 0;
%         else
%             i=i+1;  % increase time index by 1 (if n_txnode = 0)        
%         end
        
        
    end
end

%------------------------------------------------------------------


%------------------------------
% statistics
%------------------------------
n_access 
%t_access = n_access.*T_data
    % channel access time
%n_collision
%n_success
%per_user_access = n_access/(sum(n_access))

per_user_th = (n_success .* L_data * 8) / (N_slot*T_slot) / 10^6  % Mb/s
total_th = sum(per_user_th) % Mb/s7/////
fairness_index = total_th^2 / (N_user * sum(per_user_th.*per_user_th))

collision_prob = sum(n_collision)/sum(n_access)
collision_prob2 = mean(n_collision./n_access)
utilization = sum(sum(tx_state==1)) / N_slot