import matplotlib.pyplot as plt
import numpy as np

points = np.random.rand(10, 2)
fig, ax = plt.subplots()
sc = ax.scatter(points[:, 0], points[:, 1], s=100, c='blue')

dragging_idx = None  # 드래그 중인 점 인덱스

def on_press(event):
    global dragging_idx
    if event.inaxes != ax:
        return

    click = np.array([event.xdata, event.ydata])
    dists = np.linalg.norm(points - click, axis=1)
    idx = np.argmin(dists)

    if dists[idx] < 0.05:
        dragging_idx = idx

def on_motion(event):
    global dragging_idx
    if dragging_idx is None or event.inaxes != ax:
        return

    points[dragging_idx] = [event.xdata, event.ydata]
    sc.set_offsets(points)
    fig.canvas.draw_idle()

def on_release(event):
    global dragging_idx
    dragging_idx = None

# 이벤트 연결
fig.canvas.mpl_connect('button_press_event', on_press)
fig.canvas.mpl_connect('motion_notify_event', on_motion)
fig.canvas.mpl_connect('button_release_event', on_release)

plt.show()