import numpy as np
import pandas as pd
import random
import torch
import torch.nn as nn
from collections import OrderedDict

from blood import *


class Environment():
    
    def __init__(self, SETTINGS, num_bloodgroups, max_request_size, e):
        
        self.num_bloodgroups = num_bloodgroups
        self.max_request_size = max_request_size
        self.state_size = (num_bloodgroups * SETTINGS.max_age) + (num_bloodgroups * max_request_size)

        self.supply_scenario = pd.read_csv(SETTINGS.home_dir + f"supply/{SETTINGS.supply_size}/cau{round(SETTINGS.donor_eth_distr[0]*100)}_afr{round(SETTINGS.donor_eth_distr[1]*100)}_asi{round(SETTINGS.donor_eth_distr[2]*100)}_{e}.csv")
        self.demand_scenario = pd.read_csv(SETTINGS.home_dir + f"demand/{SETTINGS.avg_daily_demand}/{SETTINGS.test_days + SETTINGS.init_days}/{SETTINGS.demand_scenario}_{e}.csv")

        self.supply_index = 0

        if SETTINGS.strategy == "relimm":
            self.mismatch_weights = SETTINGS.relimm
        elif SETTINGS.strategy == "patgroups":
            self.mismatch_weights = SETTINGS.patgroup_weights

    def sample_initial_state(self, SETTINGS):
        
        inventory_vector = []

        # TODO: now we sample inventory for each day with uniform distribution, change this to sampling more young units and less old units.
        n_products = round(SETTINGS.max_inventory_size / SETTINGS.max_age)

        # Sample supply for each age upto maximum age.
        for _ in range(SETTINGS.max_age):
            inventory_vector += self.sample_supply_single_day(SETTINGS, n_products)

        demand_vector = self.sample_demand_single_day(SETTINGS)

        # Supply vector: first all bloodgroups of age 0, then all bloodgroups of age 1 etc.
        # Demand vector: first all bloodgroups with 1 unit requested, then all bloodgroups with 2 units requested etc.
        return inventory_vector + demand_vector


    def step(self, SETTINGS, state, action, day):

        if len(state) != self.state_size:
            print(f"Warning: state was expected to be of length {self.state_size}. Instead it is of length {len(state)}.")

        inventory_size = self.num_bloodgroups * SETTINGS.max_age
        inventory_vector = state[:inventory_size]
        demand_vector = state[inventory_size:]
        
        reward = self.calculate_reward(inventory_vector, demand_vector, action)
        next_state = self.next_state(SETTINGS, inventory_vector, demand_vector, action, day)

        return next_state, reward


    # TODO: calculate reward.
    def calculate_reward(self, inventory_vector, demand_vector, action):

        # print(inventory_vector)
        # print(action.tolist())
        nonexistent_products_issued = sum([min(0,i-a)*-1 for (i, a) in zip(inventory_vector, action.tolist())])
        
        return -nonexistent_products_issued


    def next_state(self, SETTINGS, inventory_vector, demand_vector, action, day):

        # Remove issued products from inventory.
        inventory_vector = [max(0,i-a) for (i, a) in zip(inventory_vector, action.tolist())]

        # Update inventory.
        outdated = inventory_vector[-self.num_bloodgroups:]                         # Get slice of inventory of age 35, which will be outdated. TODO: log outdated products somewhere.
        aging_products = inventory_vector[:-self.num_bloodgroups]                   # Get slice of inventory of age < 35.
        n_products = max(0, SETTINGS.max_inventory_size - sum(aging_products))      # Number of products to supply in order to fill complete inventory capacity.
        new_supply = self.sample_supply_single_day(SETTINGS, n_products)            # Generate new supply of age 0.
        inventory_vector = new_supply + aging_products

        # Sample new demand.
        demand_vector = self.sample_demand_single_day(SETTINGS, day)

        return inventory_vector + demand_vector


    def sample_supply_single_day(self, SETTINGS, n_products):

        # Select the next part of the supply scenario.
        data = self.supply_scenario.iloc[self.supply_index : self.supply_index + n_products]

        # Vector with a 0 for each blood group.
        vector = [0 for i in range(2**len(SETTINGS.major + SETTINGS.minor))]

        # Read each row of the selected dataframe, transform it to a bloodgroup index and add 1 to this index in the supply vector.
        for i in data.index:
            product = Blood(major = vector_to_major([data.loc[i,a] for a in SETTINGS.major]), minor = [data.loc[i,a] for a in SETTINGS.minor])
            vector[product.vector_to_bloodgroup_index()] += 1

        # Keep track of the last-read index of the supply scenario.
        self.supply_index += n_products

        return vector


    def sample_demand_single_day(self, SETTINGS, day = 0):

        # Select the part of the demand scenario belonging to the given day.
        # TODO: all requests are currently sampled as becoming known on the day of issuing. Include "Day Available" column later when desired.
        data = self.demand_scenario.loc[self.demand_scenario["Day Needed"] == day]

        demand_vector = []

        for n_units in range(self.max_request_size):
            # Vector with a 0 for each blood group.
            vector = [0 for i in range(2**len(SETTINGS.major + SETTINGS.minor))]

            # Get all rows of the data where the given number of units is requested.
            data_selected = data[data["Num Units"] == n_units + 1]

            # Read each row of the selected dataframe, transform it to a bloodgroup index and add 1 to this index in the demand vector.
            for i in data_selected.index:
                request = Blood(major = vector_to_major([data_selected.loc[i,a] for a in SETTINGS.major]), minor = [data_selected.loc[i,a] for a in SETTINGS.minor])
                vector[request.vector_to_bloodgroup_index()] += 1
            demand_vector += vector

        return demand_vector


class DQNAgent:
    
    def __init__(self, SETTINGS, num_bloodgroups, max_request_size):
        # Define DQN Layers
        self.state_vector_length = (num_bloodgroups * SETTINGS.max_age) + (num_bloodgroups * max_request_size)
        self.action_vector_length = (num_bloodgroups * SETTINGS.max_age)

        if torch.cuda.is_available():
            self.device = 'cuda' 
        # elif torch.backends.mps.is_available():
        #     self.device = 'mps'
        else:
            self.device = 'cpu'
        
        # DQN network  
        self.dqn = DQNSolver(self.state_vector_length, self.action_vector_length).to(self.device)
        self.dqn.apply(self.kaiming_init_weights)
        # self.dqn_target = DQNSolver(self.state_vector_length, self.action_vector_length).to(self.device)
        # self.dqn_target.apply(self.xavier_init_weights)
        self.optimizer = torch.optim.Adam(self.dqn.parameters(), lr=SETTINGS.alpha, weight_decay=0.0001)

        # Create memory
        # self.max_memory_size = SETTINGS.max_memory_size
        # self.STATE_MEM = torch.zeros(self.max_memory_size, self.state_vector_length)
        # self.ACTION_MEM = torch.zeros(self.max_memory_size, self.action_vector_length)
        # self.REWARD_MEM = torch.zeros(self.max_memory_size, 1)
        # self.NEXT_STATE_MEM = torch.zeros(self.max_memory_size, self.state_vector_length)
       
        # self.DONE_MEM = torch.zeros(max_memory_size, 1)
        self.ending_position = 0
        self.num_in_queue = 0
        # self.batch_size = SETTINGS.batch_size
        
        # Learning parameters
        self.gamma = SETTINGS.gamma
        self.mse = nn.SmoothL1Loss().to(self.device) # Also known as Huber loss
        self.exploration_max = SETTINGS.exploration_max
        self.exploration_rate = SETTINGS.exploration_max
        self.exploration_min = SETTINGS.exploration_min
        self.exploration_decay = SETTINGS.exploration_decay

    def dqn_load(self, SETTINGS):
        self.dqn = torch.jit.load(SETTINGS.home_dir + f"nn/{SETTINGS.model_name}.pt")
        self.dqn.eval()

    def dqn_reset(self):
        self.dqn = DQNSolver(self.state_vector_length, self.action_vector_length).to(self.device)
        self.dqn.apply(self.kaiming_init_weights);
        # self.dqn_target = DQNSolver(self.state_vector_length, self.action_vector_length).to(self.device)
        # self.dqn_target.apply(self.xavier_init_weights);

    def kaiming_init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
            m.bias.data.fill_(0.01)

    def xavier_init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            m.bias.data.fill_(0.01)

    def remember(self, state, action, reward, next_state):
        """Store the experiences in a buffer to use later"""
        self.STATE_MEM[self.ending_position] = state.float()
        self.ACTION_MEM[self.ending_position] = action.float()
        self.REWARD_MEM[self.ending_position] = reward.float()
        self.NEXT_STATE_MEM[self.ending_position] = next_state.float()
        # self.DONE_MEM[self.ending_position] = done.float()
        self.ending_position = (self.ending_position + 1) % self.max_memory_size  # FIFO tensor
        self.num_in_queue = min(self.num_in_queue + 1, self.max_memory_size)
    
    def batch_experiences(self):
        """Randomly sample 'batch size' experiences"""
        idx = random.choices(range(self.num_in_queue), k=self.batch_size)
        STATE = self.STATE_MEM[idx]
        ACTION = self.ACTION_MEM[idx]
        REWARD = self.REWARD_MEM[idx]
        NEXT_STATE = self.NEXT_STATE_MEM[idx]
        #DONE = self.DONE_MEM[idx]      
        return STATE, ACTION, REWARD, NEXT_STATE #, DONE
    
    def act(self, state):
        """Epsilon-greedy action"""
        state = torch.Tensor(state).unsqueeze(0)
        self.dqn.eval()
        with torch.no_grad():
            action_values = self.dqn(state.to(self.device))
        self.dqn.train()
        if np.random.rand() < self.exploration_rate: 
            return torch.tensor([random.randint(0,state[0][i]) for i in range(self.action_vector_length)])
        else:
            return action_values[0]
    
    def experience_replay(self):
        # Sample a batch of experiences
        STATE, ACTION, REWARD, NEXT_STATE = self.batch_experiences()
        STATE = STATE.to(self.device)
        ACTION = ACTION.to(self.device)
        REWARD = REWARD.to(self.device)
        NEXT_STATE = NEXT_STATE.to(self.device)
        # DONE = DONE.to(self.device)
        
    #     return STATE, ACTION, REWARD, NEXT_STATE

    def train(self, iters, e): 
        for _ in range(iters):    
            STATE, ACTION, REWARD, NEXT_STATE = self.experience_replay()
            self.dqn.train()
            self.dqn_target.eval()
            self.optimizer.zero_grad()

            # Q-Learning target is Q*(S, A) <- r + γ max_a Q(S', a) 
            current = self.dqn(STATE).gather(1, ACTION.long())
            with torch.no_grad():
                target = REWARD + self.gamma * self.dqn_target(NEXT_STATE).max(1).values.unsqueeze(1)
            loss = self.mse(current, target).to(self.device)
            loss.backward() # Compute gradients
            self.optimizer.step() # Backpropagate error

            # self.exploration_rate *= self.exploration_decay
            # Makes sure that exploration rate is always at least 'exploration min'
            self.exploration_rate = max(np.exp(-self.exploration_decay * e), self.exploration_min)


class DQNSolver(nn.Module):

    def __init__(self, input_shape, n_actions):
        super(DQNSolver, self).__init__()
        self.fc = nn.Sequential(OrderedDict([
          ('lin1', nn.Linear(input_shape, 8)),
          # ('sig1', nn.Sigmoid()),
          ('relu1', nn.ReLU()),
          # ('norm1', nn.BatchNorm1d(8)),
          ('lin2', nn.Linear(8, n_actions)),
          ('relu2', nn.ReLU())
        ]))    

    def forward(self, state):
        return self.fc(state)


class DQNPolicy():
    def __init__(self, dqn, device='cpu'):
        self.dqn = dqn
        self.device = device
        
    def get_action(self, state):
        state = torch.Tensor([state]).unsqueeze(0)
        self.dqn.eval()
        with torch.no_grad():
            action_values = self.dqn(state.to(self.device))
        self.dqn.train()
        action = torch.argmax(action_values).unsqueeze(0).unsqueeze(0).cpu()
        
        return action.item()

    # def get_policy(self):
    #     pi = []
    #     for s in states:
    #         a = self.get_action(s)
    #         pi.append(a)
    #     return np.array(pi)