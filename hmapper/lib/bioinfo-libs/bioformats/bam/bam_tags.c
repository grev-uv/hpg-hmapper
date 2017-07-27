#include "bam_tags.h"


//------------------------------------------------------------------------------------

bam_tag_t* bam_tag_init(const char *tag, char type, char array_type, size_t length) {
  bam_tag_t* ptr = calloc(1, sizeof(bam_tag_t));

  memset(ptr->tag, 0, 3);
  memcpy(ptr->tag, tag, 2);
  ptr->type = type;

  switch (type) {
    case BAM_TAG_TYPE_CHAR:
      ptr->data = calloc(1, sizeof(bam_char_t));
      break;
    
    case BAM_TAG_TYPE_BYTE:
      ptr->data = calloc(1, sizeof(bam_byte_t));
      break;
    
    case BAM_TAG_TYPE_UBYTE:
      ptr->data = calloc(1, sizeof(bam_ubyte_t));
      break;
    
    case BAM_TAG_TYPE_SHORT:
      ptr->data = calloc(1, sizeof(bam_short_t));
      break;
    
    case BAM_TAG_TYPE_USHORT:
      ptr->data = calloc(1, sizeof(bam_ushort_t));
      break;
    
    case BAM_TAG_TYPE_INT:
      ptr->data = calloc(1, sizeof(bam_int_t));
      break;
    
    case BAM_TAG_TYPE_UINT:
      ptr->data = calloc(1, sizeof(bam_uint_t));
      break;
    
    case BAM_TAG_TYPE_FLOAT:
      ptr->data = calloc(1, sizeof(bam_float_t));
      break;
    
    case BAM_TAG_TYPE_STRING:
      {
        ptr->data = calloc(1, sizeof(bam_string_t));
      
        bam_string_t *str = (bam_string_t*)ptr->data;
        str->length = length + 1;
        str->count = 0;
        str->data = calloc(length, sizeof(char));
      }

      break;
    
    case BAM_TAG_TYPE_HEX_STRING:
      {
        ptr->data = calloc(1, sizeof(bam_hex_string_t));
      
        bam_hex_string_t *str = (bam_hex_string_t*)ptr->data;
        str->length = length + 1;
        str->count = 0;
        str->data = calloc(length, sizeof(char));
      }

      break;

      case BAM_TAG_TYPE_BINARY:
      {
        ptr->data = calloc(1, sizeof(bam_binary_array_t));
      
        bam_binary_array_t *bin = (bam_binary_array_t*)ptr->data;
        bin->type = array_type;
        bin->length = length;
        bin->count = 0;
        
        switch (array_type) {
          case BAM_TAG_TYPE_BYTE:
            bin->data = calloc(length, sizeof(bam_byte_t));
            break;

          case BAM_TAG_TYPE_UBYTE:
            bin->data = calloc(length, sizeof(bam_ubyte_t));
            break;
          
          case BAM_TAG_TYPE_SHORT:
            bin->data = calloc(length, sizeof(bam_short_t));
            break;
          
          case BAM_TAG_TYPE_USHORT:
            bin->data = calloc(length, sizeof(bam_ushort_t));
            break;
          
          case BAM_TAG_TYPE_INT:
            bin->data = calloc(length, sizeof(bam_int_t));
            break;
          
          case BAM_TAG_TYPE_UINT:
            bin->data = calloc(length, sizeof(bam_uint_t));
            break;
          
          case BAM_TAG_TYPE_FLOAT:
            bin->data = calloc(length, sizeof(bam_float_t));
            break;
        }
      }

      break;
  }

  return ptr;
}

//------------------------------------------------------------------------------------

void bam_tag_free(bam_tag_t *ptr) {
  if (ptr) {
    if (ptr->data) {
      switch (ptr->type) {
        case BAM_TAG_TYPE_STRING:
        case BAM_TAG_TYPE_HEX_STRING:
        {
          bam_hex_string_t *str = (bam_hex_string_t*)ptr->data;
          free(str->data);
          break;
        }
        
        case BAM_TAG_TYPE_BINARY:
        {
          bam_binary_array_t *bin = (bam_binary_array_t*)ptr->data;
          free(bin->data);
          break;
        }
      }

      free(ptr->data);
    }

    free(ptr);
  }
}

//------------------------------------------------------------------------------------

void bam_tag_set_scalar(bam_tag_t *tag_p, void *data) {
  switch (tag_p->type) {
    case BAM_TAG_TYPE_BYTE:
    case BAM_TAG_TYPE_UBYTE:
      memcpy(tag_p->data, (uint8_t*)data, sizeof(uint8_t));
      break;
    
    case BAM_TAG_TYPE_SHORT:
    case BAM_TAG_TYPE_USHORT:
      memcpy(tag_p->data, (uint16_t*)data, sizeof(uint16_t));
      break;

    case BAM_TAG_TYPE_INT:
    case BAM_TAG_TYPE_UINT:
      memcpy(tag_p->data, (uint32_t*)data, sizeof(uint32_t));
      break;

    case BAM_TAG_TYPE_FLOAT:
      memcpy(tag_p->data, (float*)data, sizeof(float));
      break;
  }
}

//------------------------------------------------------------------------------------

bam_char_t bam_tag_get_char(bam_tag_t *tag_p) {
  return *(bam_char_t*)tag_p->data;
}

//------------------------------------------------------------------------------------

bam_byte_t bam_tag_get_byte(bam_tag_t *tag_p) {
  return *(bam_byte_t*)tag_p->data;
}

//------------------------------------------------------------------------------------

bam_ubyte_t bam_tag_get_ubyte(bam_tag_t *tag_p) {
  return *(bam_ubyte_t*)tag_p->data;
}

//------------------------------------------------------------------------------------

bam_int_t bam_tag_get_int(bam_tag_t *tag_p) {
  return *(bam_int_t*)tag_p->data;
}

//------------------------------------------------------------------------------------

bam_uint_t bam_tag_get_uint(bam_tag_t *tag_p) {
  return *(bam_uint_t*)tag_p->data;
}

//------------------------------------------------------------------------------------

bam_float_t bam_tag_get_float(bam_tag_t *tag_p) {
  return *(bam_float_t*)tag_p->data;
}

//------------------------------------------------------------------------------------

char* bam_tag_get_string(bam_tag_t *tag_p, size_t* length) {
  char* result = NULL;
  bam_string_t* str;

  if (tag_p->type == BAM_TAG_TYPE_STRING) {
    str = (bam_string_t*)tag_p->data;
    result = strdup(str->data);

    if (length) {
      *length = str->count + 1;
    }
  }

  return result;
}

//------------------------------------------------------------------------------------

void bam_tag_str_append(bam_tag_t *ptr, char c) {
  if (ptr->type == BAM_TAG_TYPE_STRING) {
    bam_string_t *str = (bam_string_t*)ptr->data;

    // If the character fits in the buffer, add it
    if (str->count < str->length) {
      str->data[(str->count)++] = c;
    }
  }
}

//------------------------------------------------------------------------------------

void bam_tag_str_insert(bam_tag_t *ptr, const char *str, size_t length) {
  if (ptr->type == BAM_TAG_TYPE_STRING) {
    bam_string_t *dst = (bam_string_t*)ptr->data;

    // If the string fits in the buffer, copy it
    if (dst->count + length <= dst->length) {
      memcpy(dst->data + dst->count, str, length);
      dst->count += length;
    }
  }
}

//------------------------------------------------------------------------------------

void bam_tag_str_clear(bam_tag_t *ptr, char c) {
  if (ptr->type == BAM_TAG_TYPE_STRING) {
    bam_string_t *str = (bam_string_t*)ptr->data;
    memset(str->data, c, str->length);
  }
}

//------------------------------------------------------------------------------------

char* bam_tag_to_string(const bam_tag_t **items, size_t count, size_t *out_length) {
  // Allocate a string big enough to fit all the
  // BAM tag entries in the array
  size_t result_sz = 0;
  size_t i = 0, k = 0;

  for (; i < count; ++i) {
    bam_tag_t *tag = items[i];

    // Add the size reserved for the tag name and
    // type
    result_sz += 3;

    // Add the appropriate size depending on the
    // tag type
    switch (tag->type) {
      case BAM_TAG_TYPE_CHAR:
        result_sz += sizeof(bam_char_t);
        break;
      
      case BAM_TAG_TYPE_BYTE:
        result_sz += sizeof(bam_byte_t);
        break;
      
      case BAM_TAG_TYPE_UBYTE:
        result_sz += sizeof(bam_ubyte_t);
        break;
      
      case BAM_TAG_TYPE_SHORT:
        result_sz += sizeof(bam_short_t);
        break;
      
      case BAM_TAG_TYPE_USHORT:
        result_sz += sizeof(bam_ushort_t);
        break;
      
      case BAM_TAG_TYPE_INT:
        result_sz += sizeof(bam_int_t);
        break;
      
      case BAM_TAG_TYPE_UINT:
        result_sz += sizeof(bam_uint_t);
        break;
      
      case BAM_TAG_TYPE_FLOAT:
        result_sz += sizeof(bam_float_t);
        break;
      
      case BAM_TAG_TYPE_STRING:
        {
          bam_string_t *str = (bam_string_t*)tag->data;
          result_sz += str->count + 2;
        }
        break;

      case BAM_TAG_TYPE_HEX_STRING:
      case BAM_TAG_TYPE_BINARY:
        // To-Do
        break;
    }
  }

  // Allocate the tag string
  char *result = calloc(result_sz, sizeof(char));

  // Fill the tag string data with the tag list
  for (i = 0, k = 0; i < count; ++i) {
    bam_tag_t *tag = items[i];

    // Copy the tag name
    memcpy(result + k, tag->tag, 2);
    k += 2;

    // Set the appropriate tag type and data
    switch (tag->type) {
      case BAM_TAG_TYPE_CHAR:
        {
          result[k++] = BAM_TAG_TYPE_CHAR;
          result[k++] = *((bam_char_t*)tag->data);
        }
        break;
      
      case BAM_TAG_TYPE_BYTE:
        {
          result[k++] = BAM_TAG_TYPE_BYTE;
          result[k++] = *((bam_byte_t*)tag->data);
        }
        break;
      
      case BAM_TAG_TYPE_UBYTE:
        {
          result[k++] = BAM_TAG_TYPE_UBYTE;
          result[k++] = *((bam_ubyte_t*)tag->data);
        }
        break;
      
      case BAM_TAG_TYPE_SHORT:
        {
          result[k++] = BAM_TAG_TYPE_SHORT;
          memcpy(result + k, (bam_short_t*)tag->data, sizeof(bam_short_t));
          k += sizeof(bam_short_t);
        }
        break;
      
      case BAM_TAG_TYPE_USHORT:
        {
          result[k++] = BAM_TAG_TYPE_USHORT;
          memcpy(result + k, (bam_ushort_t*)tag->data, sizeof(bam_ushort_t));
          k += sizeof(bam_ushort_t);
        }
        break;
      
      case BAM_TAG_TYPE_INT:
        {
          result[k++] = BAM_TAG_TYPE_INT;
          memcpy(result + k, (bam_int_t*)tag->data, sizeof(bam_int_t));
          k += sizeof(bam_int_t);
        }
        break;
      
      case BAM_TAG_TYPE_UINT:
        {
          result[k++] = BAM_TAG_TYPE_UINT;
          memcpy(result + k, (bam_uint_t*)tag->data, sizeof(bam_uint_t));
          k += sizeof(bam_uint_t);
        }
        break;
      
      case BAM_TAG_TYPE_FLOAT:
        {
          result[k++] = BAM_TAG_TYPE_FLOAT;
          memcpy(result + k, (bam_float_t*)tag->data, sizeof(bam_float_t));
          k += sizeof(bam_float_t);
        }
        break;
      
      case BAM_TAG_TYPE_STRING:
        {
          result[k++] = BAM_TAG_TYPE_STRING;

          bam_string_t *str = (bam_string_t*)tag->data;
          memcpy(result + k, str->data, str->count);
          result[k + str->count + 1] = '\0';
          k += str->count + 1;
        }
        break;

      case BAM_TAG_TYPE_HEX_STRING:
      case BAM_TAG_TYPE_BINARY:
        // To-Do
        break;
    }
  }

  // Return the tag string
  *out_length = result_sz;
  return result;
}

//------------------------------------------------------------------------------------

char* bam_tag_list_to_string(array_list_t *items, size_t *out_length) {
  size_t out_sz = 0;
  size_t count = array_list_size(items);
  bam_tag_t** tags = calloc(count, sizeof(bam_tag_t*));
  char* result = NULL;

  for (size_t i = 0; i < count; ++i) {
    tags[i] = array_list_get(i, items);
  }

  result = bam_tag_to_string((const bam_tag_t**)tags, count, out_length);
  free(tags);

  return result;
}

//------------------------------------------------------------------------------------

void __free_tag_cb(void* data) {
  bam_tag_t* tg = (bam_tag_t*)data;
  bam_tag_free(tg);
}

//------------------------------------------------------------------------------------

array_list_t* bam_string_to_tag(const char* bam_string, size_t length) {
  array_list_t* result = NULL;
  bam_tag_t* current = NULL;

  char tag_name[3];
  char tag_type = 'i';

  size_t i = 0, w = 0, orig_length = length;
  size_t valid_buffer = 0, is_alnum = 0;

  // Check if the first three characters of the tag string are
  // valid. If not, move the tag string until three valid characters
  // are found. If none are found, skip the tags.
  if (length > 0) {
    // Search the first character
    is_alnum = isalnum(bam_string[0]);

    if (!is_alnum) {
      while (!valid_buffer) {
        ++bam_string;
        --length;

        // Reached the end of the buffer and a valid character was not
        // found (probably the format is completely wrong)
        if (w > orig_length) {
          break;
        }

        is_alnum = isalnum(bam_string[w]);

        // If the current first character is valid, check the two following
        // characters. If all three are valid, the string is correctly formed
        // and the loop can be broken. If they are not, keep searching for a
        // valid character
        if (is_alnum) {
          if ((w + 2 < orig_length) && isalnum(bam_string[w + 1]) && isalnum(bam_string[w + 2])) {
            if (bam_string[w + 2] == BAM_TAG_TYPE_CHAR       ||
                bam_string[w + 2] == BAM_TAG_TYPE_BYTE       ||
                bam_string[w + 2] == BAM_TAG_TYPE_UBYTE      ||
                bam_string[w + 2] == BAM_TAG_TYPE_INT        ||
                bam_string[w + 2] == BAM_TAG_TYPE_UINT       ||
                bam_string[w + 2] == BAM_TAG_TYPE_FLOAT      ||
                bam_string[w + 2] == BAM_TAG_TYPE_STRING     ||
                bam_string[w + 2] == BAM_TAG_TYPE_HEX_STRING ||
                bam_string[w + 2] == BAM_TAG_TYPE_BINARY)
            valid_buffer = 1;
          }
        }
      }
    } else {
      if ((orig_length > 3) && isalnum(bam_string[1]) && isalnum(bam_string[2])) {
        valid_buffer = 1;
      }
    }
  }

  // Skip the tags if the string couldn't validated
  if (!valid_buffer) {
    return NULL;
  }

  result = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  while (i < length) {
    // Read the tag name
    memset(tag_name, 0, sizeof(tag_name));
    memcpy(tag_name, bam_string + i, 2);
    i += 2;

    // Read the tag type
    tag_type = bam_string[i++];

    // Parse the tag data using the appropriate type
    switch (tag_type) {
      case BAM_TAG_TYPE_CHAR:
        {
          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_char_t ch = bam_string[i++];
          bam_tag_set_scalar(current, &ch);
        }
        break;
      
      case BAM_TAG_TYPE_BYTE:
        {
          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_byte_t bt = bam_string[i++];
          bam_tag_set_scalar(current, &bt);
        }
        break;
      
      case BAM_TAG_TYPE_UBYTE:
        {
          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_ubyte_t ub = bam_string[i++];
          bam_tag_set_scalar(current, &ub);
        }
        break;
      
      case BAM_TAG_TYPE_SHORT:
        {
          bam_short_t sh;
          memcpy(&sh, bam_string + i, sizeof(bam_short_t));
          i += sizeof(bam_short_t);

          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_tag_set_scalar(current, &sh);
        }
        break;
      
      case BAM_TAG_TYPE_USHORT:
        {
          bam_ushort_t sh;
          memcpy(&sh, bam_string + i, sizeof(bam_ushort_t));
          i += sizeof(bam_ushort_t);

          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_tag_set_scalar(current, &sh);
        }
        break;
      
      case BAM_TAG_TYPE_INT:
        {
          bam_int_t sh;
          memcpy(&sh, bam_string + i, sizeof(bam_int_t));
          i += sizeof(bam_int_t);

          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_tag_set_scalar(current, &sh);
        }
        break;
      
      case BAM_TAG_TYPE_UINT:
        {
          bam_uint_t sh;
          memcpy(&sh, bam_string + i, sizeof(bam_uint_t));
          i += sizeof(bam_uint_t);

          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_tag_set_scalar(current, &sh);
        }
        break;
      
      case BAM_TAG_TYPE_FLOAT:
        {
          bam_float_t sh;
          memcpy(&sh, bam_string + i, sizeof(bam_float_t));
          i += sizeof(bam_float_t);

          current = bam_tag_init(tag_name, tag_type, 0, 0);
          bam_tag_set_scalar(current, &sh);
        }
        break;
      
      case BAM_TAG_TYPE_STRING:
        {
          int blen = 0;
          char* s_ptr = (char*)&bam_string[i];

          // Get the length of the BAM string
          while (*(s_ptr++)) {
            ++blen;
          }

          // Allocate the data in the BAM string
          // and fill in
          current = bam_tag_init(tag_name, tag_type, 0, blen);
          bam_tag_str_insert(current, bam_string + i, blen);
          i += blen + 1;
        }
        break;

      case BAM_TAG_TYPE_HEX_STRING:
      case BAM_TAG_TYPE_BINARY:
        // To-Do
        break;
      
      case 0x00:
        {
          // The space allocated for the tag array in the BAM file probably is
          // larger than used; otherwise the array is badly formed. If any tags
          // have been previously processed, finish the parsing. Otherwise, cancel.
          if (array_list_size(result) == 0) {
            array_list_free(result, __free_tag_cb);
            result = NULL;
          }

          return result;
        }
        break;

      default:
        {
          // The tag format is invalid. Skip the tags since the tag array is
          // probably badly formed.
          array_list_free(result, __free_tag_cb);
          return NULL;
        }
        break;
    }

    // Add the tag to the list
    array_list_insert(current, result);
  }

  return result;
}
