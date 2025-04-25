import numpy as np

class GMM:
    def __init__(self, n_components=2, max_iter=100, tol=1e-6, reg_covar=1e-6):
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.reg_covar = reg_covar

    def _gaussian_pdf(self, x, mean, cov):
        n = x.shape[1]
        cov += np.eye(n) * self.reg_covar  # regularization
        det = np.linalg.det(cov)
        inv = np.linalg.inv(cov)
        norm_const = 1.0 / (np.power(2 * np.pi, n / 2) * np.sqrt(det))
        diff = x - mean
        return norm_const * np.exp(-0.5 * np.sum(diff @ inv * diff, axis=1))

    def fit(self, X):
        n_samples, n_features = X.shape
        rng = np.random.RandomState(42)

        # 초기화
        self.weights_ = np.ones(self.n_components) / self.n_components
        self.means_ = X[rng.choice(n_samples, self.n_components, replace=False)]
        self.covariances_ = np.array([np.cov(X.T) for _ in range(self.n_components)])

        log_likelihood_old = None

        for iteration in range(self.max_iter):
            # E-step
            resp = np.zeros((n_samples, self.n_components))
            for k in range(self.n_components):
                resp[:, k] = self.weights_[k] * self._gaussian_pdf(X, self.means_[k], self.covariances_[k])
            resp_sum = np.sum(resp, axis=1, keepdims=True)
            resp /= resp_sum

            # M-step
            N_k = np.sum(resp, axis=0)
            self.weights_ = N_k / n_samples
            self.means_ = (resp.T @ X) / N_k[:, np.newaxis]

            for k in range(self.n_components):
                diff = X - self.means_[k]
                weighted_diff = resp[:, k][:, np.newaxis] * diff
                self.covariances_[k] = (weighted_diff.T @ diff) / N_k[k]
                self.covariances_[k] += np.eye(n_features) * self.reg_covar

            # 로그 우도 계산
            log_likelihood = np.sum(np.log(np.sum(resp, axis=1)))
            if log_likelihood_old is not None and abs(log_likelihood - log_likelihood_old) < self.tol:
                break
            log_likelihood_old = log_likelihood

    def predict(self, X):
        probs = np.array([
            self.weights_[k] * self._gaussian_pdf(X, self.means_[k], self.covariances_[k])
            for k in range(self.n_components)
        ])
        return np.argmax(probs, axis=0)

    def score_samples(self, X):
        probs = np.array([
            self.weights_[k] * self._gaussian_pdf(X, self.means_[k], self.covariances_[k])
            for k in range(self.n_components)
        ])
        return np.log(np.sum(probs, axis=0))
    
    
np.random.seed(0)
X1 = np.random.multivariate_normal([0, 0], np.eye(2) * 0.5, 100)
X2 = np.random.multivariate_normal([3, 3], np.eye(2) * 0.5, 100)
X = np.vstack([X1, X2])

# 학습
model = GMM(n_components=2, max_iter=1000)
model.fit(X)

# 클러스터 예측
labels = model.predict(X)
print("클러스터 레이블:", labels[:10])
import matplotlib.pyplot as plt

labelx1 = labels[:len(X1)]
labelx2 = labels[len(X1):]



plt.figure(figsize=(6, 6))
plt.scatter(X1[:, 0], X1[:, 1], c=labelx1, cmap='viridis', s=10,marker='o')
plt.scatter(X2[:, 0], X2[:, 1], c=labelx2, cmap='viridis', s=10,marker='x')
plt.scatter(model.means_[:, 0], model.means_[:, 1], c='red', s=100, marker='x')
plt.title("Custom GMM Clustering")
plt.grid(True)
plt.axis("equal")
plt.show()


from sklearn.mixture import GaussianMixture
model = GaussianMixture(2)
labels = model.fit_predict(X)
labelx1 = labels[:len(X1)]
labelx2 = labels[len(X1):]
plt.figure(figsize=(6, 6))
plt.scatter(X1[:, 0], X1[:, 1], c=labelx1, cmap='viridis', s=10,marker='o')
plt.scatter(X2[:, 0], X2[:, 1], c=labelx2, cmap='viridis', s=10,marker='x')
plt.scatter(model.means_[:, 0], model.means_[:, 1], c='red', s=100, marker='x')
plt.title("Custom GMM Clustering")
plt.grid(True)
plt.axis("equal")
plt.show()
