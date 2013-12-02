/**
 * @fileOverview Javascript for Clotho project cl02
 * @author werner@bussedesign.com
 * @version 1
 * @requires jQuery 1.9.0
 */

(function ($) {
    'use strict';
    /*jslint  plusplus: true, vars: true */
    /*global jQuery, clearTimeout, setTimeout, Image */
    
        $(function () {
        // we are using a closure to limit variable scope to this function only.
        // without this we might have a problem with firebug showing the html that is loaded after an ajax call


        function getSequence(sequenceLength) {
            var availableNumbers = [],
                randomValues = [],
                random,
                value,
                i;

            // create array of all available numbers
            for (i = 0; i < sequenceLength; i++) {
                availableNumbers[i] = i;
            }

            // in order to avoid duplicate numbers we take them out of the available array and create a new 
            // randomValues array
            for (i = 0; i < sequenceLength; i++) { // pick numbers
                // random number between 0 and the available numbers array length
                // each time we come here, availableNumbers.length will be decreased by 1 by the splice action below
                random = Math.floor(Math.random() * availableNumbers.length);
                // extract the number from the available numbers array
                // this means the length of the available numbers array will be decreased by 1
                value = availableNumbers.splice(random, 1)[0];
                // build the random numbers array
                randomValues[i] = value;
            }
            return randomValues;
        } // end getSequence





        /**
         *  Rotate a slideshow
         *  @param $slideshow - object - slideshow which needs to be updated
         *  @param newImage - object - the new image that will replace the old one
         */
        function updateSlideshow($slideshow, newImage) {
            $slideshow.prepend("<li><img src='" + newImage.src + "' /></li>");
            $slideshow.find('li:last-child').fadeOut('slow', function () {
                $(this).remove();
            });
        } // end rotateSlideshows




        /**
         *  function to update all slideshows once with a time delay between each slideshow update
         *  @param i - integer - the cyle variable. indicates where in a cycle we are
         *  @param sequence - array - with a random number sequence in which the slideshows are to be updated
         */
        function doSequence(i, sequence, imageSet) {
            var timerVar;
            
            // check for a complete sequence  
            if ((imageSet.length) === i) {
                clearTimeout(timerVar);
                return;
            }
            // get the next random number
            var nextSlideshowIndex = sequence[i],
            // update slideshow
                $thisSlideshow = $('.slideshowContainer').eq(nextSlideshowIndex),
                newImage = imageSet[i];

            updateSlideshow($thisSlideshow, newImage);
            // increment cycle variable  
            i++;
            // call the timeout function
            timerVar = setTimeout(function () {
                doSequence(i, sequence, imageSet);
            }, 500);
        } // end doSequence



        /**
         *  function to cycle all slideshows through their respective slides
         *  @param i - integer - length of cycle, e.g. how many images for all slideshows
         *  @param sequence - random sequence to access images
         *  @param numberOfImages - integer
         *  @param imageArray - array holding all image sources
         */
        function doCycle(i, sequence, numberOfImages, imageArray) {
            var timerVar2,
                imageSetIndex,
                imgSequence,
                imageSet = [];
            
            // check for a complete cycle - if cycle is done start a new one
            if (i === sequence.length) {
                clearTimeout(timerVar2);
                // start over
                i = 0;
                // create a new random sequence for all images. this is valid for one slideshow cycle
                sequence = getSequence(imageArray.length);
            }
            // generate a new random update sequence for one image set
            for (imageSetIndex = 0; imageSetIndex < numberOfImages; imageSetIndex++) {
                imageSet[imageSetIndex] = imageArray[sequence[i]];
                i++;
            }
            // update all images once
            imgSequence = getSequence(numberOfImages);
            doSequence(0, imgSequence, imageSet);
            // call the timeout function
            timerVar2 = setTimeout(function () {
                doCycle(i, sequence, numberOfImages, imageArray);
            }, 7000);
        }



        /**
         *  home page multiple slideshows
         */
        if ($('.slideshowContainer').length) {

            // create the slideshow image array
            var NUMBER_OF_IMAGES = 42, // this should be a multiple of numberOfImages
                imageArray = [],
                i = 0;

            for (i = 0; i < NUMBER_OF_IMAGES; i++) {
                imageArray[i] = new Image();
                imageArray[i].src = "assets/pictures/pic" + (i + 1) + ".jpg";
            }

            // create a random sequence for all images. this is valid for one slideshow cycle
            var sequence = getSequence(imageArray.length),
                numberOfImages = $('.slideshowContainer').length;


            // start the slideshow after initial dalay
            setTimeout(function () {
                doCycle(0, sequence, numberOfImages, imageArray);
            }, 7000);

        }

        /**
         *  Columnizer
         *  source: http://welcome.totheinter.net/columnizer-jquery-plugin/
         */
        if($('.makeColumns')) {
          $('.softwareItem').find('.description').columnize({ columns: 2, manualBreaks: true});
        }

        
        /**
         *  simple overlay
         */
        if($('#overlay').length) {
          $('.showOverlay').click(function(){
            $('#overlay').show().find('.content').show();
            return false;
          });
          
          $('#overlay').find('.close').click(function(){
            $('#overlay').hide();
          });
        }



    });
}(jQuery));