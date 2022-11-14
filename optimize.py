from RL import *

def reinforcement_learning(SETTINGS):

    mode = SETTINGS.RL_mode

    num_bloodgroups = 2**len(SETTINGS.major + SETTINGS.minor)
    max_request_size = len(SETTINGS.request_num_units_probabilities["Other"])

    agent = DQNAgent(SETTINGS, num_bloodgroups, max_request_size)
    if mode == "test":
        agent.dqn_load(SETTINGS)

    C = 0 

    for e in range(SETTINGS.episodes):
        print(f"Episode: {e}")

        env = Environment(SETTINGS, num_bloodgroups, max_request_size, e)
        state = env.sample_initial_state(SETTINGS)

        C += 1
        for day in range(1, SETTINGS.init_days + 1)
            if day % 50 == 0:
                print(f"Day {day}")

                action = agent.act(state)
                state_next, reward = env.step(SETTINGS, state, torch.tensor([round(a) for a in action.tolist()]), day)
                # agent.remember(torch.Tensor([state]), action, torch.tensor([reward]).unsqueeze(0), torch.Tensor([state_next]))
                state = state_next

        for day in range(SETTINGS.init_days + 1, SETTINGS.init_days + SETTINGS.test_days + 1):
            if day % 50 == 0:
                print(f"Day {day}")

            action = agent.act(state)
            state_next, reward = env.step(SETTINGS, state, torch.tensor([round(a) for a in action.tolist()]), day)
            agent.remember(torch.Tensor([state]), action, torch.tensor([reward]).unsqueeze(0), torch.Tensor([state_next]))
            state = state_next

        if mode == "train":
            mean_rewards.append(mean_reward)
            agent.train(train_iters, e)
        
            if C % SETTINGS.nn_update_iter == 0:
                print(f'NN update {C}')
                agent.dqn_target = agent.dqn
         
            if e % 5 == 0:
                model_scripted = torch.jit.script(agent.dqn) # Export to TorchScript
                model_scripted.save(SETTINGS.home_dir + "nn/test.pt") # Save
