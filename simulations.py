import numpy as np

class SnakesAndLadders():
    '''
    Simulate the Snakes and Ladders game.

    Parameters
    ----------
    snakes : list of tuples
        snakes indicates the snakes, e.g. [(bottom, top), (bottom, top), ...]
    ladders : list of tuples
        ladders indicates the ladders, e.g. [(bottom, top), (bottom, top), ...]
    '''
    def __init__(self, snakes: list, ladders: list):
        self.P = self._create_P_matrix(snakes, ladders)
        np.random.seed(42)

    def _create_P_matrix(self, snakes, ladders):
        P = np.zeros((101, 101))
        for i in range(101):
            prob = 1/6
            P[i, (i+1):(min((i+6), 100) + 1)] = prob
            P[i, i] += 1 - P[i, :].sum()
        for snake in snakes:
            end, start = snake
            P[:, end] += P[:, start]
            P[:, start] = 0
        for ladder in ladders:
            start, end = ladder
            P[:, end] += P[:, start]
            P[:, start] = 0
        return P
    
    def simulation(self):
        '''Simulate the entire game 1000 times and
        return the mean and std of number of dice throws.'''
        n = []
        for _ in range(1000):
            n_round = self._round()
            n.append(n_round)
        return self._report(n)

    def _round(self):
        X = 0
        n = 0
        while X != 100:
            X = self._dice(self.P[X, :])
            n += 1
            if n >= 10000:
                raise ValueError('Too many iterations withou reaching the end. Review the snakes.')
        return n

    def _dice(self, p):
        result = np.random.multinomial(1, p)
        return np.argmax(result)

    def _report(self, n):
        return [np.mean(n), np.std(n)]


class MMInfinity():
    '''
    Simulate the row M/M/infinity.

    Parameters
    ----------
    lambda : float
        The rate at which people arrive
    mu : float
        The rate of time service
    '''
    def __init__(self, lambda_, mu):
        self.lambda_ = lambda_
        self.mu = mu
        np.random.seed(42)

    def simulation(self):
        '''
        Simulate 1000 times the row.
        
        Returns
        -------
        time_mean : float
            The mean of time the row takes to become empty again.
        time_std : float
            The std. of time the row takes to become empty again.
        max_mean : float
            The mean of the maximum size achieved by the row.
        max_std : float
            The std. of the maximum size achieved by the row.
        '''
        times = []
        max_sizes = []
        for _ in range(1000):
            time, max_size = self._round()
            times.append(time)
            max_sizes.append(max_size)

        return self._report(times, max_sizes)

    def _round(self):
        X = 1
        max_X = X
        time = self._rexp(self.lambda_)
        while X != 0:
            forward = self._rexp(self.lambda_)
            backward = self._rexp(X*self.mu)
            time += min(forward, backward)
            if backward < forward:
                X -= 1
            else:
                X += 1
            if X > max_X:
                max_X = X
        return time, max_X

    def _report(self, times, max_sizes):
        time_mean = np.mean(times)
        time_std = np.std(times)
        max_mean = np.mean(max_sizes)
        max_std = np.std(max_sizes)
        return time_mean, time_std, max_mean, max_std

    def _rexp(self, rate):
        u = np.random.uniform(0, 1)
        rexp = -np.log(1 - u)/rate
        return rexp
        # return np.random.exponential(1/rate)

class Brownian():
    '''
    Simulate the Brownian motion.

    Parameters
    ----------
    sigma_2 : float
        The diffusion coefficient
    drift : float
        The drift
    start : float
        Start of Brownian motion
    end : float
        End of Brownian motion
    delta_t : float
        The partition size of [start, end]
    M : float
        Value to be calculated the median of time to it be achieved
    '''
    def __init__(self, sigma_2=1, drift=0, start=0, end=10, delta_t=0.01, M=-1.5):
        self.sigma_2 = sigma_2
        self.drift = drift
        self.start = start
        self.end = end
        self.delta_t = delta_t
        self.M = M
        np.random.seed(42)

    def simulation(self):
        '''
        Simulate 1000 times the Brownian motion.
        
        Returns
        -------
        max_90_quantile : float
            The 0.9-quantile of maximum of Brownian motion
        tau_median : float
            The median of time to achieve M
        '''
        max_Bs = []
        taus = []
        for _ in range(1000):
            max_B, tau = self._round()
            max_Bs.append(max_B)
            taus.append(tau)

        return self._report(max_Bs, taus)

    def _round(self):
        B = 0
        max_B = B
        tau = 0
        ts = np.arange(self.start + self.delta_t, self.end + self.delta_t, self.delta_t)
        B_old = B
        M_positive = self.M > 0
        for t in ts:
            normal = np.random.normal()
            B = B_old + np.sqrt(self.delta_t)*np.sqrt(self.sigma_2)*normal + self.delta_t*self.drift
            if B > max_B:
                max_B = B
            if M_positive == True and tau == 0 and B > self.M:
                tau = t
            if M_positive == False and tau == 0 and B < self.M:
                tau = t
            B_old = B
        if tau == 0:
            tau = 10
        return max_B, tau

    def _report(self, max_Bs, taus):
        max_90_quantile = np.quantile(max_Bs, 0.9)
        tau_median = np.median(taus)
        return max_90_quantile, tau_median


def main():
    print('Exercise 1')
    print()

    print('Exercise a)')
    a = SnakesAndLadders([], [])
    a_mean, a_std = a.simulation()
    print('Mean = {:6.3f}, std = {:6.3f}'.format(a_mean, a_std))
    print()

    print('Exercise b)')
    b = SnakesAndLadders(
        [(6, 16), (11, 49), (19, 62), (24, 87), (26, 47), (60, 64), (73, 93), (75, 95), (78, 98)],
        [(1, 38), (4, 14), (9, 31), (21, 42), (28, 84), (36, 44), (51, 67), (71, 91), (80, 100)])
    b_mean, b_std = b.simulation()
    print('Mean = {:6.3f}, std = {:6.3f}'.format(b_mean, b_std))
    print()

    print('Exercise c)')
    c = SnakesAndLadders(
        [(26, 47), (78, 98)],
        [(1, 38), (4, 14), (9, 31), (21, 42), (28, 84), (36, 44), (51, 67), (71, 91), (80, 100)])
    c_mean, c_std = c.simulation()
    print('Mean = {:6.3f}, std = {:6.3f}'.format(c_mean, c_std))
    print()

    print('Exercise 2')
    print()

    print('Exercise a)')
    row = MMInfinity(lambda_=1.4, mu=1.3)
    time_mean, time_std, max_mean, max_std = row.simulation()
    print('Mean = {:6.3f}, std = {:6.3f}'.format(time_mean, time_std))
    print()

    print('Exercise b)')
    print('Mean = {:6.3f}, std = {:6.3f}'.format(max_mean, max_std))
    print()

    print('Exercise 3')
    print()

    print('Exercise a)')
    no_drift = Brownian()
    max_90_quantile, tau_median = no_drift.simulation()
    print('Max value 90-quantile = {:6.3f}'.format(max_90_quantile))
    print()

    print('Exercise c)')
    print('Median of time to achieve M = {:6.3f}'.format(tau_median))
    print()

    print('Exercise e)')
    drift = Brownian(drift=-0.14)
    max_90_quantile, tau_median = drift.simulation()
    print('Max value 90-quantile = {:6.3f}'.format(max_90_quantile))
    print()

    print('Exercise f)')
    print('Median of time to achieve M = {:6.3f}'.format(tau_median))
    print()

    
if __name__ == '__main__':
    main()
