# A-Hybrid-Method-for-Pavement-Crack-Width-Measurement
This repository contains the hybrid method for pavement crack width measurement. The details can be found here: https://www.sciencedirect.com/science/article/abs/pii/S026322412200505X#:~:text=The%20hybrid%20method%20obtains%20the,close%20to%20the%20orthogonal%20direction.

Measurement samples of our method <br>
Original image       |  Binary image    |  Measurement
:-------------------------:|:-------------------------:|:-------------------------:
![33](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/4af87d33-466d-4a62-a9cb-a7a621667f5f) | ![33](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/0784fb34-24db-41e5-85fc-1683b24050e7) | ![image](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/6b7f7c14-5d8f-4a25-92d0-fc130bf15156)
![338](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/6f7e5170-a04f-4abc-8169-9c043b3c927c) | ![338](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/4855e1cf-06af-4314-878b-2fb7eb79d31e) | ![image](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/5c58773b-1557-4efe-b378-5dcbb5ed24df)
![34](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/203c606f-5147-42ff-9fb3-5e1b5dcbfda2) |![34](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/c07f8ac5-e602-408c-a975-e19bbfa1a7f7) | ![image](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/06e8f857-413b-49a9-bbbd-0e5f0a6ca860)
![31](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/41a8349a-3773-4603-859d-0acb8efe8ec6)| ![31](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/48c1ff35-bb1f-406a-90a8-67bad06ca537) |![image](https://github.com/JeremyOng96/A-Hybrid-Method-for-Pavement-Crack-Width-Measurement/assets/17587452/bbd418f9-6aa7-482e-8eec-d27722b160fb)


The benefits of the hybrid method are as follows:
   <ul>
     <li>Allows the user to tune the strength of the orthogonal projection and the shortest method</li>
     <li>Prevents overestimation effects by the orthogonal projection method when cracks are curved</li>
     <li>More distinct crack width measurements than the shortest method</li>
     <li>More robust to binary images that have unparallel boundaries</li>
   </ul>
   
How to use:
   <ol>
    <li>Input binary image</li>
    <li>Tune the pruning strength</li>
    <li>Tune the balancing coefficient</li>
   </ol>


   
