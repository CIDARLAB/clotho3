'use strict';

/**
 * @name unique
 * @example {{items | unique:attribute | json}}
 * Filters out all duplicate items from an array by checking the specified key
 * @param [key] {string} the name of the attribute of each object to compare for uniqueness
 if the key is empty, the entire object will be compared
 if the key === false then no filtering will be performed
 * @return {array}
 */
Application.Interface.filter('unique', ['$parse', function ($parse) {

    return function (items, filterOn) {

        if (filterOn === false) {
            return items;
        }

        if ((filterOn || angular.isUndefined(filterOn)) && angular.isArray(items)) {
            var hashCheck = {}, newItems = [],
                get = angular.isString(filterOn) ? $parse(filterOn) : function (item) { return item; };

            var extractValueToCompare = function (item) {
                return angular.isObject(item) ? get(item) : item;
            };

            angular.forEach(items, function (item) {
                var valueToCheck, isDuplicate = false;

                for (var i = 0; i < newItems.length; i++) {
                    if (angular.equals(extractValueToCompare(newItems[i]), extractValueToCompare(item))) {
                        isDuplicate = true;
                        break;
                    }
                }
                if (!isDuplicate) {
                    newItems.push(item);
                }

            });
            items = newItems;
        }
        return items;
    };
}]);

/**
 * @name highlight
 * @example

 <pre class="prettyprint">
 &lt;label&gt;&lt;input type=&quot;checkbox&quot; ng-model=&quot;caseSensitive&quot;&gt; Case Sensitive?&lt;/label&gt;
 &lt;input placeholder=&quot;Enter some text to highlight&quot; value=&quot;you&quot; ng-model=&quot;highlightText&quot;&gt;
 &lt;p ng-bind-html-unsafe=&quot;&#x27;Hello there, how are you today? I\'m fine thank you.&#x27; | highlight:highlightText:caseSensitive&quot;&gt;&lt;/p&gt;

 &lt;style&gt;
 .ui-match { background: yellow; }
 &lt;/style&gt;
 </pre>

 * @param text {string} haystack to search through
 * @param search {string} needle to search for
 * @param [caseSensitive] {boolean} optional boolean to use case-sensitive searching
 */
Application.Interface.filter('highlight', function () {
    return function (text, search, caseSensitive) {
        if (search || angular.isNumber(search)) {
            text = text.toString();
            search = search.toString();
            if (caseSensitive) {
                return text.split(search).join('<span class="ui-match">' + search + '</span>');
            } else {
                return text.replace(new RegExp(search, 'gi'), '<span class="ui-match">$&</span>');
            }
        } else {
            return text;
        }
    };
});