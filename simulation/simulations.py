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
    
if __name__ == '__main__':
    main()
