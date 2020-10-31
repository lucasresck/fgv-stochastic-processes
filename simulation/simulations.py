import numpy as np

class SnakesAndLadders():
    '''
    Simulate the Snakes and Ladders game.

    Parameters
    ----------
    snakes : list of tuples
        snakes indicates the snakes, e.g. [(start, end), (start, end), ...]
    ladders : list of tuples
        ladders indicates the ladders, e.g. [(start, end), (start, end), ...]
    '''
    def __init__(self, snakes: list, ladders: list):
        self.P = self._create_P_matrix(snakes, ladders)

    def _create_P_matrix(self, snakes, ladders):
        P = np.zeros((101, 101))
        for i in range(101):
            prob = 1/6
            P[i, (i+1):(min((i+6), 100) + 1)] = prob
            P[i, i] += 1 - P[i, :].sum()
        return P
    
    def simulation(self):
        '''Simulate the entire game 1000 times and
        return the mean and variance of number of dice throws.'''
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
        return n

    def _dice(self, p):
        result = np.random.multinomial(1, p)
        return np.argmax(result)

    def _report(self, n):
        return [np.mean(n), np.var(n)]
    
