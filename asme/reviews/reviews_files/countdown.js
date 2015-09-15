/**************************************************
countdown.js
Huibert Kwakernaak
February 2007
**************************************************/

function countDown(){
   var Countdown = showcountdown(countdown);
   box.value = Countdown;
   if (countdown < 5*60 && warn == 'no'){
      box.style.backgroundColor = '#FFFF00';
      box.style.color = 'FF0000';
      warn = 'yes';
   }
   if (countdown > 0){
      setTimeout(countDown,1000);
      countdown = countdown-1;
   }
   else{
      box.style.backgroundColor = '#FF0000';
      box.style.color = 'FFFFFF';
   }
}

function showcountdown(countdown){
   var hours = Math.floor(countdown/3600);
   countdown = countdown-hours*3600;
   var minutes = Math.floor(countdown/60);
   var seconds = countdown-minutes*60;
   minutes = minutes.toString();
   seconds = seconds.toString();
   if (minutes.length < 2){minutes = '0' + minutes;}
   if (seconds.length < 2){seconds = '0' + seconds;}
   var Countdown = minutes + ':' + seconds;
   if (hours > 0){Countdown = hours + ':' + Countdown;}
   return Countdown;
}

