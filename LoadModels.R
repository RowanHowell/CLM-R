createNetwork("Activity0.txt", "Localization0.txt", "FEARnet0.txt", "LocSpecificActivity0.txt", "model0", directory = "networkData")
model0 = loadNetwork("model0OE.txt")
ICs0 = makeICs("MetaAIC0.csv", "MetaLIC0.csv", model0, directory = "networkData")
createNetwork("Activity0.txt", "Localization0Cdc157A.txt", "FEARnet0.txt", "LocSpecificActivity0.txt", "model0Cdc157A", directory = "networkData")
model0Cdc157A = loadNetwork("model0Cdc157AOE.txt")
createNetwork("Activity0.txt", "Localization0Mob12A.txt", "FEARnet0.txt", "LocSpecificActivity0.txt", "model0Mob12A", directory = "networkData")
model0Mob12A = loadNetwork("model0Mob12AOE.txt")
createNetwork("Activity0.txt", "Localization0DPM.txt", "FEARnet0.txt", "LocSpecificActivity0.txt", "model0DPM", directory = "networkData")
model0DPM = loadNetwork("model0DPMOE.txt")

createNetwork("Activity1.txt", "Localization1.txt", "FEARnet0.txt", "LocSpecificActivity1.txt", "model1", directory = "networkData")
model1 = loadNetwork("model1OE.txt")
ICs1 = makeICs("MetaAIC1.csv", "MetaLIC1.csv", model1, directory = "networkData")
createNetwork("Activity1.txt", "Localization1Cdc157A.txt", "FEARnet0.txt", "LocSpecificActivity1Cdc157A.txt", "model1Cdc157A", directory = "networkData")
model1Cdc157A = loadNetwork("model1Cdc157AOE.txt")
createNetwork("Activity1.txt", "Localization1Mob12A.txt", "FEARnet0.txt", "LocSpecificActivity1.txt", "model1Mob12A", directory = "networkData")
model1Mob12A = loadNetwork("model1Mob12AOE.txt")
createNetwork("Activity1.txt", "Localization1DPM.txt", "FEARnet0.txt", "LocSpecificActivity1Cdc157A.txt", "model1DPM", directory = "networkData")
model1DPM = loadNetwork("model1DPMOE.txt")

createNetwork("Activity2.txt", "Localization2.txt", "FEARnet0.txt", "LocSpecificActivity2.txt", "model2", directory = "networkData")
model2 = loadNetwork("model2OE.txt")
ICs2 = makeICs("MetaAIC2.csv", "MetaLIC2.csv", model2, directory = "networkData")
createNetwork("Activity2.txt", "Localization2Cdc157A.txt", "FEARnet0.txt", "LocSpecificActivity2Cdc157A.txt", "model2Cdc157A", directory = "networkData")
model2Cdc157A = loadNetwork("model2Cdc157AOE.txt")
createNetwork("Activity2.txt", "Localization2Mob12A.txt", "FEARnet0.txt", "LocSpecificActivity2.txt", "model2Mob12A", directory = "networkData")
model2Mob12A = loadNetwork("model2Mob12AOE.txt")
createNetwork("Activity2.txt", "Localization2DPM.txt", "FEARnet0.txt", "LocSpecificActivity2Cdc157A.txt", "model2DPM", directory = "networkData")
model2DPM = loadNetwork("model2DPMOE.txt")


createNetwork("Activity3.txt", "Localization3.txt", "FEARnet0.txt", "LocSpecificActivity3.txt", "model3", directory = "networkData")
model3 = loadNetwork("model3OE.txt")
ICs3 = makeICs("MetaAIC3.csv", "MetaLIC3.csv", model3, directory = "networkData")
createNetwork("Activity3.txt", "Localization3DPM.txt", "FEARnet0.txt", "LocSpecificActivity3Cdc157A.txt", "model3DPM", directory = "networkData")
model3DPM = loadNetwork("model3DPMOE.txt")

createNetwork("Activity3a.txt", "Localization3a.txt", "FEARnet0.txt", "LocSpecificActivity3a.txt", "model3a", directory = "networkData")
model3a = loadNetwork("model3aOE.txt")
ICs3a = makeICs("MetaAIC5.csv", "MetaLIC5.csv", model3a, directory = "networkData")

createNetwork("Activity3a.txt", "Localization3aDPM.txt", "FEARnet0.txt", "LocSpecificActivity3aCdc157A.txt", "model3aDPM", directory = "networkData")
model3aDPM = loadNetwork("model3aDPMOE.txt")

createNetwork("Activity4.txt", "Localization4.txt", "FEARnet0.txt", "LocSpecificActivity4.txt", "model4", directory = "networkData")
model4 = loadNetwork("model4OE.txt")
ICs4 = makeICs("MetaAIC3.csv", "MetaLIC3.csv", model4, directory = "networkData")

createNetwork("Activity4a.txt", "Localization4a.txt", "FEARnet0.txt", "LocSpecificActivity4a.txt", "model4a", directory = "networkData")
model4a = loadNetwork("model4aOE.txt")
ICs4a = makeICs("MetaAIC3.csv", "MetaLIC3.csv", model4a, directory = "networkData")

createNetwork("Activity5.txt", "Localization5.txt", "FEARnet0.txt", "LocSpecificActivity5.txt", "model5", directory = "networkData")
model5 = loadNetwork("model5OE.txt")
ICs5 = makeICs("MetaAIC5.csv", "MetaLIC5.csv", model5, directory = "networkData")

createNetwork("Activity5.txt", "Localization5DPM.txt", "FEARnet0.txt", "LocSpecificActivity5Cdc157A.txt", "model5DPM", directory = "networkData")
model5DPM = loadNetwork("model5DPMOE.txt")
ICs5DPM = makeICs("MetaAIC5.csv", "MetaLIC5.csv", model5DPM, directory = "networkData")

createNetwork("Activity5.txt", "Localization5Cdc157A.txt", "FEARnet0.txt", "LocSpecificActivity5Cdc157A.txt", "model5Cdc157A", directory = "networkData")
model5Cdc157A = loadNetwork("model5Cdc157AOE.txt")

createNetwork("Activity5.txt", "Localization5a.txt", "FEARnet0.txt", "LocSpecificActivity5.txt", "model5a", directory = "networkData")
model5a = loadNetwork("model5aOE.txt")

createNetwork("Activity5.txt", "Localization5aDPM.txt", "FEARnet0.txt", "LocSpecificActivity5Cdc157A.txt", "model5aDPM", directory = "networkData")
model5aDPM = loadNetwork("model5aDPMOE.txt")
ICs5a = makeICs("MetaAIC5.csv", "MetaLIC5.csv", model5aDPM, directory = "networkData")



