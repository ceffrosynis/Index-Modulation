# SI-MM-OFDM-IM

## Introduction

Index modulation (IM) refers to a family of modulation techniques that rely on the activation states of some resources/building blocks for information embedding. The resources/building blocks can be either physical, e.g., antenna, subcarrier, time slot, and frequency carrier, or virtual, e.g., virtual parallel channels, signal constellation, space-time matrix, and antenna activation order. A distinct feature of IM is that part of the information is implicitly embedded into the transmitted signal.

## What is in this repo?

In this repo we implement three of the basics IM schemes:
* **Single Mode OFDM Index Modulation (SM-OFDM-IM)**
* **Dual Mode OFDM Index Modulation (DM-OFDM-IM)**
* **Multiple Mode OFDM Index Modulation (MM-OFDM-IM)**

and we propose and implement a novel (IM) scheme, with a much more improved Bit Error Rate (BER) performance, the
* **Sub-block Index Multiple Mode OFDM Index Modulation (SI-MM-OFDM-IM)**.

![alt text](https://github.com/ceffrosynis/Index-Modulation/blob/master/images/Arrow%20Diagram%20Casual%20Strcture(2).png)

This repo also contains a recearch paper
* **Obfuscation_of_Index_Modulated_OFDM_Signals.pdf**

which contains 
* **Model Architecture** 
* **BER Performance Evaluation (Unencrypted/Encrypted)**
* **Modulation Obfuscation Security Evaluation (Unencrypted/Encrypted)**


## Implementation

* **Transmitter**

* **Channel Model**
  * **Rician Fading Channel**
  * **Perfect Channel Estimation**

* **Receiver**
  * **Detector**
    * **Viterbi-like Algorithm**
    * **Better performance than the Optimal ML Detector**


## SM-OFDM-IM

## DM-OFDM-IM

## MM-OFDM-IM

## SI-OFDM-IM
