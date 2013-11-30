Application.Foundation.service('LocalNgStorage', ['$window', 'PubSub', function($window, PubSub) {

    //defaults
    var refStorage = $window.localStorage;
    var serializer = JSON;
    //prefix to prevent name-clashes
    var prefix = "clotho_";


    // --- Check for support ---
    // FUTURE - add this in later to be more robust
    // e.g. https://github.com/grevory/angular-local-storage/blob/master/localStorageModule.js
    // and could fallback to cookies for small strings (<4kb)

    // Checks the browser to see if local storage is supported
    var browserSupportsLocalStorage = function () {
        try {
            return ('localStorage' in window && window['localStorage'] !== null);
        } catch (e) {
            return false;
        }
    };

    // Checks the browser to see if cookies are supported
    // look into angular $cookies if we want to implement this as a backup
    // note the limitations of cookies before committing to that
    var browserSupportsCookies = function() {
        try {
            return navigator.cookieEnabled ||
                ("cookie" in document && (document.cookie.length > 0 ||
                    (document.cookie = "test").indexOf.call(document.cookie, "test") > -1));
        } catch (e) {
            return false;
        }
    };

    // ---- functionality ---

    var clear = function() {
        refStorage.clear();
        return( this );
    };

    // returns an item, or optional defaultValue if not found
    var getItem = function( key, defaultValue ){
        var value = refStorage.getItem( prefix+key );
        if (value == null){
            // check to see if default is falsy
            return((typeof( defaultValue ) != "undefined") ? defaultValue : null );
        } else {
            return(serializer.parse( value ) );

        }
    };

    // return a boolean whether a given object exists
    var hasItem = function( key ){
        return( refStorage.getItem( prefix+key ) != null );
    };

    //remove a given item from storage
    var removeItem = function( key ){
        refStorage.removeItem( prefix+key );
        return( this );
    };

    //adds an item to storage, automatically serializing.
    //NOTE -  cannot add functions or private variables
    var setItem = function( key, value ){
        //prevent adding of empty objects
        if (!value && value!==0 && value!=="") return false;
        refStorage.setItem(prefix+key, serializer.stringify( value ) );
        return( this );
    };

    // ----- LOCAL STORAGE EVENTS ----
    /* NOTES
     this half is for model update broadcasting across pages (via localStorage)
     the other half is using PUB/SUB (within pages + to bubble up updates)

     we can do something like $scope.$watch($window('storage'), function() {} )

     note - implementation of localStorage events varies & is unreliable
     - event handlers only invoked for current window / tab where data is written / deleted (even though spec states all windows)
     - may not be called when a key is updated, and only when it is created / deleted
     - storage deletion only returns key, not deleted value
     */

    var handle_storage_change = function(e) {
        //console.log("change made to local storage");
        if (!e) { e = window.event; }

        //TODO - better checking for e.key
        var uuid = e.key.replace(prefix, '') || '';
        console.log("handle_storage_event");
        PubSub.trigger("model_change", uuid);
        PubSub.trigger("model_change:"+uuid, getItem(uuid));
    };

    if (window.addEventListener) {
        window.addEventListener("storage", handle_storage_change, false);
        console.log("local storage Listener added");
    } else {
        window.attachEvent("onstorage", handle_storage_change); //IE8
        console.log("local storage Listener added for IE8");
    }

    // ------- FACADE -----------
    return {
        getPrefix : function() {return prefix},
        isSupported : browserSupportsLocalStorage,
        clear : clear,
        getItem : getItem,
        hasItem : hasItem,
        removeItem : removeItem,
        setItem : setItem
    }

}]);