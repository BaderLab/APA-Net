from setuptools import setup, find_packages

setup(
    name='APA-Net',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'torch',
        'scipy',
        'wandb',
        'tqdm',
    ],
    entry_points={
        'console_scripts': [
            'apa-train=apamodel.train_script:main',
        ],
    },
)
