from functools import partial

import tqdm


def tqdm_callback(N, notebook=True):
    """
    Return a :module:`tqdm` progress bar expecting *N* iterations,
    either suitable with jupyter if *notebook* is true and for the
    terminal otherwise. The progress bar includes an additional method
    :function:`callback` (function of one ignored parameter) meant to
    be past as a callback function called to update the progress bar.
    """
    if notebook:
        progress_bar = tqdm.tqdm_notebook(total=N)
    else:
        progress_bar = tqdm.tqdm(total=N)
    def callback(self, i):
        self.update()
    progress_bar.callback = partial(callback, progress_bar)
    return progress_bar


if __name__ == '__main__':
    from time import sleep

    N = 10

    def test(N, callback):
        for i in range(N):
            sleep(0.2)
            callback(i)

    with tqdm_callback(N, notebook=False) as progress_bar:
        test(N, progress_bar.callback)
